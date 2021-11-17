#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mu(double xx, double yy){
  NumericVector out(2);
  out[0] = -xx; out[1] = -yy;
  return(out);}


// // [[Rcpp::export]]
// IntegerVector ranking(NumericVector x){
//   Function f("rank");
//   return f(x);
// }

// [[Rcpp::export]]
IntegerVector ranking(NumericVector invec){
  int leng = invec.size();
  IntegerVector y(leng);
  for(int i=0; i<leng; ++i){
    y[i] = sum(invec<invec[i]) + 1;
  }
  return(y);
}


// // [[Rcpp::export]]
// NumericVector order_cpp(NumericVector invec){
//   int leng = invec.size();
//   NumericVector y = clone(invec);
//   for(int i=0; i<leng; ++i){
//     y[sum(invec<invec[i])] = i+1;
//   }
//   return(y);
// }
// // [[Rcpp::export]]
// NumericVector ranking(NumericVector invec){
//   return order_cpp(order_cpp(invec));
// }


// [[Rcpp::export]]
NumericVector segcirc(NumericVector p1, NumericVector p2,double xc, double yc, double r){
  double dx = p2[0]-p1[0];
  double dy = p2[1]-p1[1];
  double m = dy/dx;
  double n = p1[1]-m*p1[0];
  double A = 1+pow(m,2);
  double B = 2*(m*n-m*yc-xc);
  double C = pow(xc,2)+pow(n,2)-2*n*yc+pow(yc,2)-pow(r,2);
  //# Intersection point (x,y)
  //# http://apetrilla.blogspot.com.uy/2011/12/circle-and-line-intersection.html
  double x = (-B+sqrt(pow(B,2)-4*A*C))/(2*A);
  double y = m*x+n;
  //# Intersection point (xb,yb)
  double xb = (-B-sqrt(pow(B,2)-4*A*C))/(2*A);
  double yb = m*xb+n;
  NumericVector salida_segcirc;
  //# Determine whether (x,y) or (xb,yb) lies on the segment [p1,p2]
  if(p1[0]!=p2[0]){
    NumericVector c = {p1[0],p2[0],x};
    IntegerVector ranki = ranking(c);
    // IntegerVector ranki = match(c,c.sort());
    NumericVector c2 = {p1[0],p2[0],xb};
    // IntegerVector ranki2 = match(c2,c2.sort());
    IntegerVector ranki2 = ranking(c2);
    if(ranki[2]==2){
      salida_segcirc = {x,y};
    }
    else if(ranki2[2]==2){
      salida_segcirc = {xb,yb};
    }
  //return salida_segcirc;
  }
  else{
    NumericVector c3 = {p1[1],p2[1],y};
    // IntegerVector ranki3 = match(c3,c3.sort());
    IntegerVector ranki3 = ranking(c3);
    NumericVector c4 = {p1[1],p2[1],yb};
    // IntegerVector ranki4 = match(c4,c4.sort());
    IntegerVector ranki4 = ranking(c4);
    if(ranki3[2]==2){
      salida_segcirc = {x,y};
    }
    else if(ranki4[2]==2){
      salida_segcirc = {xb,yb};
    }
  }
  return salida_segcirc;
  }

// [[Rcpp::export]]
NumericVector segelli(NumericVector p1, NumericVector p2,double a, double b){
  double dx = p2[0]-p1[0];
  double dy = p2[1]-p1[1];
  double m = dy/dx;
  double dr = sqrt(pow(b,2)*pow(dx,2)+pow(a,2)*pow(dy,2));
  double D = p2[0]*p1[1] - p1[0]*p2[1];
  // # Intersection point (x,y)
  double x = ( -pow(a,2)*m*D/dx + a*b*sqrt(pow(dr/dx,2)-pow(D/dx,2)))/(pow(dr/dx,2)) ;
  double y = ( pow(b,2)*D/dx + a*b*dy/dx*sqrt(pow(dr/dx,2)-pow(D/dx,2)))/(pow(dr/dx,2)) ;
  // # Intersection point (xb,yb)
  double xb = ( -pow(a,2)*m*D/dx - a*b*sqrt(pow(dr/dx,2)-pow(D/dx,2)))/(pow(dr/dx,2)) ;
  double yb = ( pow(b,2)*D/dx - a*b*dy/dx*sqrt(pow(dr/dx,2)-pow(D/dx,2)))/(pow(dr/dx,2)) ;
  // # Determine whether (x,y) or (xb,yb) lies on the segment [p1,p2]
  NumericVector salida_segelli(2);
  if(p1[0]!=p2[0]){
    NumericVector cb = {p1[0],p2[0],x};
    IntegerVector rankib = ranking(cb);
    NumericVector cb2 = {p1[0],p2[0],xb};
    IntegerVector rankib2 = ranking(cb2);
    if(rankib[2]==2){
      salida_segelli = {x,y};
    } else if(rankib2[2]==2){
      salida_segelli = {xb,yb};
    }
  } else{
    NumericVector cb3 = {p1[1],p2[1],y};
    IntegerVector rankib3 = ranking(cb3);
    NumericVector cb4 = {p1[1],p2[1],yb};
    IntegerVector rankib4 = ranking(cb4);
    if(rankib3[2]==2){
      salida_segelli = {x,y};
    } else if(rankib4[2]==2){
      salida_segelli = {xb,yb};
    }
  }
  return salida_segelli;}


// // [[Rcpp::export]]
// NumericVector segelli(NumericVector p1, NumericVector p2,double a, double b){
//   double dx = p2[0]-p1[0];
//   double dy = p2[1]-p1[1];
//   double m = dy/dx;
//   double dr = sqrt(pow(b,2)*pow(dx,2)+pow(a,2)*pow(dy,2));
//   double D = p2[0]*p1[1] - p1[0]*p2[1];
//   // # Intersection point (x,y)
//   double x = ( -pow(a,2)*m*D/dx + a*b*sqrt(pow(dr/dx,2)-pow(D/dx,2)))/(pow(dr/dx,2)) ;
//   double y = ( pow(b,2)*D/dx + a*b*dy/dx*sqrt(pow(dr/dx,2)-pow(D/dx,2)))/(pow(dr/dx,2)) ;
//   // # Intersection point (xb,yb)
//   double xb = ( -pow(a,2)*m*D/dx - a*b*sqrt(pow(dr/dx,2)-pow(D/dx,2)))/(pow(dr/dx,2)) ;
//   double yb = ( pow(b,2)*D/dx - a*b*dy/dx*sqrt(pow(dr/dx,2)-pow(D/dx,2)))/(pow(dr/dx,2)) ;
//   // # Determine whether (x,y) or (xb,yb) lies on the segment [p1,p2]
//   NumericVector salida_segelli(2);
//   if(p1[0]!=p2[0]){
//     NumericVector cb = {p1[0],p2[0],x};
//     IntegerVector rankib = match(cb,cb.sort());
//     NumericVector cb2 = {p1[0],p2[0],xb};
//     IntegerVector rankib2 = match(cb2,cb2.sort());
//     if(rankib[2]==2){
//       salida_segelli = {x,y};
//     } else if(rankib2[2]==2){
//       salida_segelli = {xb,yb};
//     }
//   } else{
//     NumericVector cb3 = {p1[1],p2[1],y};
//     IntegerVector rankib3 = match(cb3,cb3.sort());
//     NumericVector cb4 = {p1[1],p2[1],yb};
//     IntegerVector rankib4 = match(cb4,cb4.sort());
//     if(rankib3[2]==2){
//       salida_segelli = {x,y};
//     } else if(rankib4[2]==2){
//       salida_segelli = {xb,yb};
//     }
//   }
// return salida_segelli;}
    

// [[Rcpp::export]]
NumericMatrix rbmd_cpp(
    int N, double h,double a, double b, double r, double xc,
    double yc, double sigma,NumericVector ini, double eps){
  
  NumericVector x(N);
  NumericVector y(N);
  
  NumericVector xr;
  NumericVector yr;
  
  // NumericVector norm1(N); norm1 = Rcpp::rnorm(N,0,sigma*sqrt(h));
  // NumericVector norm2(N); norm2 = Rcpp::rnorm(N,0,sigma*sqrt(h));
  
  //NumericVector xa(N);
  // std::partial_sum(norm1.begin(), norm1.end(), norm1.begin());
  
  x[0] = ini[0]; y[0] = ini[1];
  xr[0] = x[0]; yr[0] = y[0];
  
  for (int i = 0; i < N; i++) {
    double xaux=x[i]+sigma*sqrt(h)*R::rnorm(0,1)+h*mu(x[i],y[i])[0];
    double yaux=y[i]+sigma*sqrt(h)*R::rnorm(0,1)+h*mu(x[i],y[i])[1];
    
    if ( pow(xaux-xc,2)+pow(yaux-yc,2) < pow(r,2) ) {
      // # 1) si esta dentro del círculo
      // # rebotar fuera del círculo y dentro de la elipse
      
      NumericVector inte = segcirc({x[i],y[i]},{xaux,yaux},xc,yc,r);
      double di = sqrt(pow(xaux-inte[0],2)+pow(yaux-inte[1],2)) ;
      
      NumericVector u = {x[i]-inte[0],y[i]-inte[1]};
      u = u/sqrt(sum(pow(u,2))) ;
      NumericVector v = {xc-inte[0],yc-inte[1]};
      v = v/sqrt(sum(pow(v,2))) ;
      double ang = acos(sum(u*v)) ;
      
      NumericVector vrot(2);
      vrot[0] =  cos(ang) * v[0] + sin(ang) * v[1];
      vrot[1] = -sin(ang) * v[0] + cos(ang) * v[1];
    
    if (sum(pow(vrot-u,2))<eps){
      vrot[0] = cos(ang) * v[0] - sin(ang) * v[1];
      vrot[1] = sin(ang) * v[0] + cos(ang) * v[1];
    }
    
    double xaux1 = inte[0]+di*vrot[0] ;
    double yaux1 = inte[1]+di*vrot[1] ;
    
      // # si me da un punto dentro del circulo o fuera de la elipse que lo tire
      // # sino que lo guarde
    if((pow(xaux1-xc,2)+pow(yaux1-yc,2) < pow(r,2)) | (pow(b,2)*pow(xaux1,2)+pow(a,2)*pow(yaux1,2) > pow(a,2)*pow(b,2))){ 
      x[i+1] = x[i];
      y[i+1] = y[i];
      // is.intc<-c(is.intc,0) # cuenta los rebotes para fuera del circulo
      // is.intintc<-c(is.intintc,1) # cuenta los desechos por rebotar dentro del circulo o fuera de la elipse
      // is.inte<-c(is.inte,0)  # cuenta los rebotes para dentro de la elipse
      // is.intinte<-c(is.intinte,0) # cuenta los desechos por rebotar fuera de la elipse o fuera del círculo
      // is.noint<-c(is.noint,0)
    } else{
      x[i+1] = xaux1;
      y[i+1] = yaux1;
      // Rcout << "A ver 2 " << i << std::endl;
      xr.push_back(inte[0]);xr.push_back(x[i+1]);
      yr.push_back(inte[1]);yr.push_back(y[i+1]);
      // xr.push_back(inte[0]);xr.push_back(xaux1);
      // yr.push_back(inte[1]);yr.push_back(yaux1);
      // is.intc<-c(is.intc,1) # cuenta los rebotes para fuera del circulo
      // is.intintc<-c(is.intintc,0) # cuenta los desechos por rebotar dentro del circulo o fuera de la elipse
      // is.inte<-c(is.inte,0)  # cuenta los rebotes para dentro de la elipse
      // is.intinte<-c(is.intinte,0) # cuenta los desechos por rebotar fuera de la elipse o fuera del círculo
      // is.noint<-c(is.noint,0)
    }
    
    } 
    else if (pow(b,2)*pow(xaux,2)+pow(a,2)*pow(yaux,2) > pow(a,2)*pow(b,2)) { //# Change in center !=0
    // # 2) si esta fuera de la elipse
    // # rebotar dentro de la elipse y fuera del círculo
    
    NumericVector inte = segelli({x[i],y[i]},{xaux,yaux},a,b);
    double di = sqrt(pow(xaux-inte[0],2)+pow(yaux-inte[1],2));
    
    NumericVector u = {x[i]-inte[0],y[i]-inte[1]};
    u = u/sqrt(sum(pow(u,2))) ;
    NumericVector v = {pow(b,2)*inte[0],pow(a,2)*inte[1]};
    v = v/sqrt(sum(pow(v,2))) ;
    double ang = acos(sum(u*v)) ;
    
    NumericVector vrot(2);
    vrot[0] =  cos(ang) * v[0] + sin(ang) * v[1];
    vrot[1] = -sin(ang) * v[0] + cos(ang) * v[1];
    
    if (sum(pow(vrot-u,2))<eps){
      vrot[0] = cos(ang) * v[0] - sin(ang) * v[1];
      vrot[1] = sin(ang) * v[0] + cos(ang) * v[1];
    }
    
    double xaux2 = inte[0]+di*vrot[0];
    double yaux2 = inte[1]+di*vrot[1];
    
    // #si me da un punto fuera de la elipse o dentro del círculo que lo tire
    // # sino que lo guarde
    if ((pow(b,2)*pow(xaux2,2)+pow(a,2)*pow(yaux2,2) > pow(a,2)*pow(b,2) ) |( pow(xaux2-xc,2)+pow(yaux2-yc,2) < pow(r,2))){
      x[i+1] = x[i];
      y[i+1] = y[i];
      // is.intc<-c(is.intc,0) # cuenta los rebotes para fuera del circulo
      // is.intintc<-c(is.intintc,0) # cuenta los desechos por rebotar dentro del circulo o fuera de la elipse
      // is.inte<-c(is.inte,0)  # cuenta los rebotes para dentro de la elipse
      // is.intinte<-c(is.intinte,1) # cuenta los desechos por rebotar fuera de la elipse o fuera del círculo
      // is.noint<-c(is.noint,0)
    } else{
      x[i+1] = xaux2;
      y[i+1] = yaux2;
      xr.push_back(inte[0]);xr.push_back(x[i+1]);
      yr.push_back(inte[1]);yr.push_back(y[i+1]);
      // xr.push_back(inte[0]);xr.push_back(xaux2);
      // yr.push_back(inte[1]);yr.push_back(yaux2);
      // is.intc<-c(is.intc,0) # cuenta los rebotes para fuera del circulo
      // is.intintc<-c(is.intintc,0) # cuenta los desechos por rebotar dentro del circulo o fuera de la elipse
      // is.inte<-c(is.inte,1)  # cuenta los rebotes para dentro de la elipse
      // is.intinte<-c(is.intinte,0) # cuenta los desechos por rebotar fuera de la elipse o fuera del círculo
      // is.noint<-c(is.noint,0)
    }
    
    }
    else{
    // # 3) Guardar el punto porque está dentro de la elipse y fuera del círculo
    x[i+1] = xaux;
    y[i+1] = yaux;
    xr.push_back(x[i+1]);
    yr.push_back(y[i+1]);
    // is.intc<-c(is.intc,0) # cuenta los rebotes para fuera del circulo
    // is.intintc<-c(is.intintc,0) # cuenta los desechos por rebotar dentro del circulo o fuera de la elipse
    // is.inte<-c(is.inte,0)  # cuenta los rebotes para dentro de la elipse
    // is.intinte<-c(is.intinte,0) # cuenta los desechos por rebotar fuera de la elipse o fuera del círculo
    // is.noint<-c(is.noint,1)
  }
      
  }
  NumericMatrix salida = cbind(xr,yr);
  return salida;  
}


