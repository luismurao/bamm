// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp ;

// Function to create the list of species present
// at each site (row) of the PAM.
// @param m The presence absence matrix

Rcpp::List spList(Rcpp::NumericMatrix m) {

  int R = m.nrow();
  int C = m.cols();
  Rcpp::List hp(R);
  Rcpp::IntegerVector vec =  Rcpp::seq(1, C);
  for(int i =0; i< R; i++){
    Rcpp::NumericVector a = m.row(i);
    hp[i] = vec[a>0];
  }
  return hp;
}

// Function to concatetante (combine) two integer vectors
// @param a Integer vector
// @param b Integer vector


Rcpp::IntegerVector conc(IntegerVector a, IntegerVector b){
  std::vector<int> a_b = as<std::vector<int> >(a);
  a_b.insert(a_b.end(),b.begin(),b.end());
  IntegerVector a_b_conc = wrap(a_b);
  return a_b_conc;
}

// Function to sample a integer vector
// @param x Integer vector
// @param size Number of elements to sample
// @param replace Logical. Sample with replacement

IntegerVector csample_num( IntegerVector x,
                           int size,
                           bool replace,
                           NumericVector prob = NumericVector::create()
) {
  IntegerVector ret = RcppArmadillo::sample(x, size, replace, prob);
  return ret;
}

// Function to find the intersection between two vectors.
// @param x Integer vector
// @param y Integer vector

Rcpp::IntegerVector intersectx(IntegerVector x, IntegerVector y){

  std::vector<int> v2 = as<std::vector<int> >(x);
  std::vector<int> v1 = as<std::vector<int> >(y);
  std::sort(v1.begin(), v1.end());
  std::sort(v2.begin(), v2.end());
  std::vector<int> v_intersection;

  std::set_intersection(v1.begin(), v1.end(),
                        v2.begin(), v2.end(),
                        std::back_inserter(v_intersection));

  IntegerVector inter = wrap(v_intersection);
  return inter;
}
// Function to find the difference between two vectors.
// @param x Integer vector
// @param y Integer vector

IntegerVector std_setdiff(IntegerVector x, IntegerVector y) {

  std::vector<int> v2 = as<std::vector<int> >(x);
  std::vector<int> v1 = as<std::vector<int> >(y);
  std::vector<int> out;

  std::set_difference(x.begin(), x.end(), y.begin(), y.end(),
                      std::inserter(out, out.end()));
  IntegerVector sdiff = wrap(out);
  return sdiff;
}


// Function to create the list of randomized entries (species) for
// each row in the PAM.
// @param m Presence-Absence-Matrix (PAM) or a binary matrix with columns
//          representing species and rows sites.
// @param niter Number of itereations to permute the PAM.
// @return Returns a list of length equal to the number of sites
// in the PAM (n=nrows(m)). Each entrie of the list has the permuted
// species.

// [[Rcpp::export]]
Rcpp::List rList(Rcpp::NumericMatrix m, int niter) {

  int R = m.nrow();
  //int C = m.cols();
  Rcpp::List hp = spList(m);

  Rcpp::IntegerVector hpI =  Rcpp::seq(0, R-1);
  for (int i =0; i< niter; i++){
    //Rcpp::IntegerVector AB = RcppArmadillo::sample<IntegerVector>(hpI, 2, false);
    Rcpp::IntegerVector AB = csample_num(hpI, 2, false);
    int arand=AB[0];
    int brand = AB[1];

    IntegerVector a =  hp[arand];
    IntegerVector b = hp[brand];
    IntegerVector ab = intersectx(a,b) ;
    int l_ab=ab.size();
    int l_a=a.size();
    int l_b=b.size();
    //int inab = l_ab == l_a || l_ab == l_b;
    if ((l_ab == l_a || l_ab == l_b)==false){
      IntegerVector a_b = conc(a,b);
      IntegerVector tot0 = setdiff(a_b,ab);
      int l_tot0=tot0.size();
      Rcpp::IntegerVector tot1 = csample_num(tot0, l_tot0, false);
      int L=l_a-l_ab;

      Rcpp::IntegerVector inHp0 = conc(ab,tot1[Rcpp::seq(0, L-1)]);
      Rcpp::IntegerVector inHp1=conc(ab,tot1[Rcpp::seq(L, l_tot0-1)]);
      hp[AB[0]] = inHp0;
      hp[AB[1]] = inHp1;
    }


  }

  return hp;
}


Rcpp::NumericMatrix fill_matrix(Rcpp::NumericMatrix m, Rcpp::List sps){
  int R = m.nrow();
  int C = m.cols();
  Rcpp::NumericMatrix m1(R,C);
  for(int i =0;i<R;i++){
    //Rcpp::NumericVector reng  = m.row(i);
    Rcpp::NumericVector spsin = sps[i];
    for(int k =0; k<spsin.size();k++){
      m1(i,spsin[k]-1) =1;

    }
  }
  return m1;
}


// Function to permute a PAM.
// @param m Presence-Absence-Matrix (PAM) or a binary matrix with columns
//          representing species and rows sites.
// @param niter Number of itereations to permute the PAM.
// @return Returns a permuted PAM.
// [[Rcpp::export]]

Rcpp::NumericMatrix permute_matrix(Rcpp::NumericMatrix m, int niter){
  List sps = rList(m,niter);
  Rcpp::NumericMatrix m1 = fill_matrix(m,sps);
  return m1;
}


Rcpp::NumericVector Quantile(Rcpp::NumericVector x, Rcpp::NumericVector probs) {
  // implementation of type 7
  const size_t n=x.size(), np=probs.size();
  if (n==0) return x;
  if (np==0) return probs;
  Rcpp::NumericVector index = (n-1.)*probs, y=x.sort(), x_hi(np), qs(np);
  Rcpp::NumericVector lo = Rcpp::floor(index), hi = Rcpp::ceiling(index);

  for (size_t i=0; i<np; ++i) {
    qs[i] = y[lo[i]];
    x_hi[i] = y[hi[i]];
    if ((index[i]>lo[i]) && (x_hi[i] != qs[i])) {
      double h;
      h = index[i]-lo[i];
      qs[i] = (1.-h)*qs[i] + h*x_hi[i];
    }
  }
  return qs;
}



// Function to permute a categorize the null dispersion distribution

// @param dfield A matrix with entries representing the observed values of dispersion field.
// @param dfield_rand A matrix with columns representing random dispersion field distribution
// @param lower_interval Lower value of the distrbution
// @param upper_interval Upper value of the distrbution
// @resturn A categorical vector
// [[Rcpp::export(null_dispersion_field_cat)]]

Rcpp::NumericVector null_dispersion_field_cat(Rcpp::NumericMatrix dfield,
                                              Rcpp::NumericMatrix dfield_rand,
                                              double lower_interval,
                                              double upper_interval){
  //qqq=Quantile(pam_vals[h,],prob=1:20/20);
  int nrows = dfield_rand.nrow();
  //int ncols = dfield_rand.ncol();
  Rcpp::NumericVector dfield_cat(nrows);
  NumericVector probs = {lower_interval,upper_interval};
  Rcpp::NumericVector q_distfield(2);
  for(int i=0; i< nrows; i++){
    //Rcpp::NumericVector dfield_i = dfield_rand.row(i);
    q_distfield = Quantile(dfield_rand.row(i),probs);
    double dfield_val = dfield.row(i)[0];
    if(dfield_val <= q_distfield[0]){
      dfield_cat[i] =1;
    }
    if(dfield_val >= q_distfield[1]){
      dfield_cat[i] =3;
    }

  }
  return dfield_cat;
}
