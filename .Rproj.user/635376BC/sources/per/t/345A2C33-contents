// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <Rmath.h>
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

//IntegerVector std_setdiff(IntegerVector x, IntegerVector y) {

//  std::vector<int> v2 = as<std::vector<int> >(x);
//  std::vector<int> v1 = as<std::vector<int> >(y);
//  std::vector<int> out;

//  std::set_difference(x.begin(), x.end(), y.begin(), y.end(),
//                      std::inserter(out, out.end()));
//  IntegerVector sdiff = wrap(out);
//  return sdiff;
//}


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
// Fastball algorithm
// Function to create the list of randomized entries (species) for
// each row in the PAM.
// @param inputList presence-absence-list created using the function spList.
//        The list contains the species present in each site of the PAM.
// @param numSwaps Number of iterations to permute the PAM.
// @return Returns a list of length equal to the number of sites
// in the PAM (n=nrows(m)). Each entree of the list has the permuted species
// @details The "fastball" algorithm is described in Godard and Neal (2022) and
//          this is implementation is taken from
//          \url{https://github.com/zpneal/fastball/blob/main/fastball.cpp}.
//          Please when using the function for publications cite
//          Godard K, Neal ZP (2022). "fastball: a fast algorithm to randomly
//          sample bipartite graphs with fixed degree sequences."
//          Journal of Complex Networks, 10(6), cnac049.

Rcpp::List fastball_cpp(Rcpp::List inputList, int numSwaps) {

  //get number of rows
  int numRows = inputList.length();

  //convert input list into a 2D std::vector
  std::vector<std::vector<int>> oneLocs (numRows);
  for(int i = 0; i < numRows; i++) {
    oneLocs[i] = Rcpp::as<std::vector<int> > (inputList[i]);
  }

  //conduct row swaps a total of (numSwaps) times
  for (int i = 1; i <= numSwaps; i++) {

    //get two random row numbers
    int r1Index = R::runif(0,1) * numRows;
    int r2Index = r1Index;
    while (r2Index == r1Index) {
      r2Index = R::runif(0,1) * numRows;
    }

    //create references to the two rows being mixed
    std::vector<int> & r1 = oneLocs[r1Index];
    std::vector<int> & r2 = oneLocs[r2Index];
    if (r1.size() == 0 || r2.size() == 0) {
      continue;
    }

    //generate iterators for first pass through rows
    std::vector<int>::iterator first1 = r1.begin();
    std::vector<int>::iterator last1 = r1.end();
    std::vector<int>::iterator first2 = r2.begin();
    std::vector<int>::iterator last2 = r2.end();
    size_t intersectionLength = 0;

    //find the length of the intersection
    while (first1!=last1 && first2!=last2)
    {
      if (*first1<*first2) ++first1;
      else if (*first2<*first1) ++first2;
      else {
        intersectionLength += 1;
        ++first1; ++first2;
      }
    }

    //calculate length of symmetric difference
    size_t r1SymDiffSize = r1.size() - intersectionLength;
    size_t r2SymDiffSize = r2.size() - intersectionLength;
    size_t symDiffSize = r1SymDiffSize + r2SymDiffSize;

    if (symDiffSize == 0) {
      continue;
    }

    //create vector of zeros and ones
    //represents which row elements of the symmetric difference are placed in
    std::vector<int> swapLocations (symDiffSize);
    std::fill(swapLocations.begin(), swapLocations.begin() + r1SymDiffSize, 0);
    std::fill(swapLocations.begin() + r1SymDiffSize, swapLocations.end(), 1);

    //shuffle swapLocations using Fisher-Yates shuffle
    for (size_t i = 0; i < swapLocations.size() - 1; i++) {
      size_t j = i + R::runif(0,1) * (swapLocations.size() - i);
      std::swap(swapLocations[i],swapLocations[j]);
    }

    //create vectors to store output of curveball swaps
    std::vector<std::vector<int> > curveballRows (2);
    curveballRows[0].reserve(r1.size());
    curveballRows[1].reserve(r2.size());

    //generate iterators for sweep through r1 and r2
    first1 = r1.begin();
    last1 = r1.end();
    first2 = r2.begin();
    last2 = r2.end();
    std::vector<int>::iterator swapIterator = swapLocations.begin();

    //compare elements in r1 and r2 until end of a vector is reached
    while (first1!=last1 && first2!=last2)
    {
      //element in row1 is less than row2
      if (*first1<*first2) {
        //use swapLocations to add element to n1 or n2
        curveballRows[*swapIterator].push_back(*first1);
        //increment iterators
        ++swapIterator;
        ++first1;
      }
      //element in row1 is greater than row2
      else if (*first2<*first1) {
        //use swapLocations to add element to n1 or n2
        curveballRows[*swapIterator].push_back(*first2);
        //increment iterators
        ++swapIterator;
        ++first2;
      }
      //element in row1 is equal to row2
      else {
        //add element to poth arrays and increment both iterators
        curveballRows[0].push_back(*first1);
        curveballRows[1].push_back(*first2);
        ++first1; ++first2;
      }
    }

    //pass through remainder of r1
    while (first1 != last1) {
      curveballRows[*swapIterator].push_back(*first1);
      ++swapIterator;
      ++first1;
    }

    //pass through remainder of r2
    while (first2 != last2) {
      curveballRows[*swapIterator].push_back(*first2);
      ++swapIterator;
      ++first2;
    }

    //clear the rows for r1 and r2 in oneLocs
    r1.clear();
    r2.clear();

    //create references to the shuffled rows
    std::vector<int> & newV1 = curveballRows[0];
    std::vector<int> & newV2 = curveballRows[1];

    //insert the data for the shuffled rows back into oneLocs
    r1.insert(r1.end(), newV1.begin(), newV1.end());
    r2.insert(r2.end(), newV2.begin(), newV2.end());
  }

  //Return randomized adjacency list
  Rcpp::List randomizedList (numRows);

  for (int i = 0; i < numRows; i++) {
    Rcpp::IntegerVector temp = Rcpp::wrap(oneLocs[i]);
    randomizedList[i] = temp;
  }
  return randomizedList;
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

// Function to permute a PAM.
// @param m Presence-Absence-Matrix (PAM) or a binary matrix with columns
//          representing species and rows sites.
// @param niter Number of itereations to permute the PAM.
// @return Returns a permuted PAM.
// [[Rcpp::export]]

Rcpp::NumericMatrix permute_matrix_fb(Rcpp::NumericMatrix m, int niter){
  Rcpp::List hp = spList(m);
  List sps = fastball_cpp(hp,niter);
  Rcpp::NumericMatrix m1 = fill_matrix(m,sps);
  return m1;
}


/*This takes an integer n and returns one random integer
  This function originally if from the package picante

*/
int intrand(int n)
{
  double u;
  u = unif_rand();
  return((int)(u*n));/*Cast the double as an integer*/
}

// Function to permute a PAM using the independent swap algorithm
// which is implemented in the picante package. Please if you
// use this function cite
// S.W. Kembel, P.D. Cowan, M.R. Helmus, W.K. Cornwell, H. Morlon,
// D.D. Ackerly, S.P. Blomberg, and C.O. Webb. 2010. Picante: R tools
// for integrating phylogenies and ecology. Bioinformatics 26:1463-1464.

// @param m Presence-Absence-Matrix (PAM) or a binary matrix with columns
//          representing species and rows sites.
// @param niter Number of iterations to permute the PAM.
// @return Returns a permuted PAM.

// [[Rcpp::export]]
Rcpp::NumericMatrix permute_matrix_indswap(Rcpp::NumericMatrix pam, int niter){
  int swap;
  int swapped;
  int i,j,k,l;
  int row = pam.nrow();
  int column = pam.ncol();
  double tmp;
  for(swap= 0; swap < niter; swap++){
    swapped = 0;
    while(swapped == 0){
      i = intrand(row);
      while((j=intrand(row))==i);
      k =  intrand(column);
      while((l = intrand(column)) == k);
      if((pam(i,k)>0 && pam(j,l)>0 && pam(i,l)+pam(j,k)==0)||(pam(i,k)+pam(j,l)==0 && pam(i,l)>0 && pam(j,k)>0))
      {
        tmp = pam(i,k);
        pam(i,k) = pam(j,k);
        pam(j,k) = tmp;
        tmp = pam(i,l);
        pam(i,l) = pam(j,l);
        pam(j,l) = tmp;
        swapped = 1;
      }
    }
  }
  return pam;
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

