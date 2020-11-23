#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]


arma::vec dloss(arma::vec y, arma::mat x, double gamma, arma::vec beta, arma::vec xbj, arma::mat xj){
   int n = y.size();
   int l = xj.n_cols;
   arma::vec ei(n), wi(n), pi(n), dlossv(l);
   arma::vec eta = x*beta + xbj;
   arma::mat xjt = xj.t();
   for(int i=0; i<n; i++){
	   ei(i) = exp((1+gamma)*eta(i));
	   wi(i) = std::pow(std::pow(ei(i),y(i))/(1+ei(i)),gamma/(gamma+1));
	   pi(i) = gamma*wi(i)*( y(i) - ei(i)/(1+ei(i)) )/n;
   }
   dlossv = - xjt*pi;
   return dlossv;
}



double objloss(arma::vec y, arma::mat x, double gamma, arma::vec beta, arma::vec xbj){
   int n = y.size();
   arma::vec eta = x*beta + xbj;
   arma::vec ei(n), wi(n), pi(n);
   double objlossv = 0;
   for(int i=0; i<n; i++){
	   ei(i) = exp((1+gamma)*eta(i));
	   wi(i) = std::pow(std::pow(ei(i),y(i))/(1+ei(i)),gamma/(gamma+1));
	   objlossv = objlossv - wi(i)/n;
   }
   return objlossv;
}





double MCP(double theta, double MCPlam, double MCPa){
  double thetaabs = std::abs(theta);
  if(thetaabs > MCPlam*MCPa){
    return std::pow(MCPa,2)*MCPlam/2;
  }else {
    return MCPlam*thetaabs - std::pow(thetaabs,2)/(2*MCPa);
  }
}


double objMCP(arma::vec theta, double MCPlam, double MCPa){
  int l = theta.size();
  arma::vec Ptheta(l);
  double objlossv = 0;
  for(int i = 0; i<l; i++){
	  Ptheta(i) = MCP(std::abs(theta(i)),MCPlam,MCPa);
	  objlossv = objlossv + Ptheta(i);
  }
  return objlossv;
}


double groupsingle(arma::vec b1, double MCPlamgroup, double MCPlamsingle, double MCPa){
   int l = b1.size();
   double penalty = MCP(norm(b1,2), MCPlamgroup, MCPa) + objMCP(b1.rows(1,l-1), MCPlamsingle, MCPa);
   return penalty;
}


double groupsinglevector(arma::vec b1, double MCPlamgroup, double MCPlamsingle, double MCPa, int p){
   int l = b1.size();
   int q = l/(p+1)-1;
   double groupsinglep = 0;
   for (int i=0; i<p; i++){
	   groupsinglep = groupsinglep + groupsingle(b1.rows((q+1)*i,(q+1)*i+q), MCPlamgroup, MCPlamsingle, MCPa);
   }
   return groupsinglep;
}


double dMCP(double theta, double MCPlam, double MCPa){
  double thetaabs = std::abs(theta);
  if(thetaabs > MCPlam*MCPa){
    return 0;
  }else {
    return MCPlam - thetaabs/MCPa;
  }
}



double sgn(double x){
  double norm = std::abs(x);
  if(x < 0){
    return -x/norm;
  }else if(x > 0){
    return x/norm;
  }else{
	return 0;
  }

}


arma::vec dMCPvector(arma::vec theta, double MCPlam, double MCPa){
  int l = theta.size();
  arma::vec Ptheta(l);
  for(int i = 0; i<l; i++){
	  Ptheta(i) = dMCP(std::abs(theta(i)),MCPlam,MCPa);
  }
  return Ptheta;
}


arma::vec amj(arma::vec y, arma::mat x, double gamma, double val0, arma::vec b0, arma::vec HG,
         double MCPlamgroup, double MCPlamsingle, double MCPa, arma::vec xbj, int pp){
	int l = b0.size();
	int amji = 0;
	int amj_num = 20;
	arma::vec amjresults(l+2);
	arma::vec b1;
	double val1;

    while(amji < amj_num){
		b1 = b0 - std::pow(0.5,amji) * HG;
        val1 = objloss(y, x, gamma, b1, xbj) + groupsingle(b1, MCPlamgroup, MCPlamsingle, MCPa);

		amji = amji + 1;
		if (val1 < val0){
			break;
		}
	}
	amjresults(0) = val1;
	amjresults(1) = amji;
	amjresults.rows(2,l+1) = b1;
	return amjresults;
}



// [[Rcpp::export]]

List gammalogisticsgMCP (arma::vec y, arma::mat x, arma::mat w, double gamma, arma::vec bhat0,
         double MCPlamgroup, double MCPlamsingle, double MCPa, int max_iter, double eps, int sear){
	int n = y.size();
    int q = x.n_cols;
    int p = w.n_cols/(q+1);
	int iter = 0;
	int l = bhat0.size();
	int max_iter2 = 20;
	arma::vec iter2(p);
	iter2.zeros();
	double d = 1;
	double dd = 1;
	arma::mat I1 = arma::ones<arma::mat>(n,1);
	x.insert_cols(0,I1);
	arma::mat xj0(n,q+1), xj(n,q+1);
	arma::mat XW = arma::join_rows(x, w);
	arma::vec b0(l), b1(l), xb(n), b1gi(q+1), b1j(q+1), xbj(n), HG(q+1), b00(q+1), bhat(l), amj_iter(p+1);
	arma::vec amjresults;
	double val1;
	double norm_bg, dp_group;
	b1 = bhat0;
	while(iter < max_iter && d > eps && max(b1) < 3){
		iter = iter + 1;
		b0 = b1;
		xb = XW * b1;

		xj = XW.cols(0,q);
		b1j = b1.rows(0,q);
		xbj = xb - xj * b1j + xj0 * b1gi;
		val1 = objloss(y, xj, gamma, b1j, xbj) + groupsingle(b1j, sqrt(q+1)*MCPlamgroup, MCPlamsingle, MCPa);
		b00 = b1.rows(0,q);
		HG = dloss(y, xj, gamma, b00, xbj, xj);
		if (sear == 0) {
            // Armijo search determining the step size
            amjresults = amj(y, xj, gamma, val1, b00, HG, sqrt(q+1)*MCPlamgroup, MCPlamsingle, MCPa, xbj, 0);
            val1 = amjresults(0);
            amj_iter(0) = amjresults(1);
            b1gi = amjresults.rows(2,q+2);
        }
		if (sear == 1) {
            // Given an appropriate step size
            b1gi = b00 - 0.001*HG;
		}
		b1.rows(0,q) = b1gi;
		xj0 = XW.cols(0,q);

		for (int j=1; j<p+1; j++){
			xj = XW.cols((q+1)*j,(q+1)*j+q);
			b1j = b1.rows((q+1)*j,(q+1)*j+q);
			xbj = xb - xj * b1j + xj0 * b1gi;
			val1 = objloss(y, xj, gamma, b1j, xbj) + groupsingle(b1j, sqrt(q+1)*MCPlamgroup, MCPlamsingle, MCPa);

			while(iter2(j-1) < max_iter2){
				b00 = b1j;
				iter2(j-1) = iter2(j-1) + 1;
				HG = dloss(y, xj, gamma, b00, xbj, xj);
				norm_bg = std::max(norm(b00,2),0.0001);
				dp_group = dMCP(norm_bg, sqrt(q+1)*MCPlamgroup, MCPa);
				HG(0) = HG(0) + dp_group * b00(0) / norm_bg;
				for (int k=1; k<q+1; k++){
					HG(k) = HG(k) + dp_group * b00(k) / norm_bg + dMCP(b00(k), MCPlamsingle, MCPa) * sgn(b00(k));
				}

				amjresults = amj(y, xj, gamma, val1, b00, HG, sqrt(q+1)*MCPlamgroup, MCPlamsingle, MCPa, xbj, 0);
				val1 = amjresults(0);amj_iter(j) = amjresults(1);
				b1gi = amjresults.rows(2,q+2);

				b1.rows((q+1)*j,(q+1)*j+q) = b1gi;
				dd=norm(b1gi-b00,2);
				if (dd < eps/10 or max(b1gi) > 3){
					break;
				}
				b1j = b1gi;
		    }
			xj0 = XW.cols((q+1)*j,(q+1)*j+q);
		}
		d=norm(b1.rows(q+1,(q+1)*p+q)-b0.rows(q+1,(q+1)*p+q),2);
	}
	xbj.zeros();
	for(int i=0; i<l; i++){
		if (std::abs(b1(i)) < 0.001){
			b1(i) = 0;
		}
	}
	bhat = b1;
	double objlossvalue = objloss(y, XW, gamma, bhat, xbj);
	List res = List::create(Named("beta")=bhat.rows(1,l-1),
                          _["beta0"]= bhat.row(0), _["obj"] = objlossvalue, _["iter"] = iter, _["iter2"] = iter2);
    return(res);
}

