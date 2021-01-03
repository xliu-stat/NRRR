// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List pinv_my(arma::mat Bg){
    /* the same as ginv() in R*/
    arma::mat U;
    arma::vec s;
    arma::mat V;
    bool pinv_success = true;
    int p = Bg.n_cols;
    int q = Bg.n_rows;
    arma::mat inv_mat(p, q);
    int i;
    int l, r=0;
    double tol;

    //bool svd_success = svd(U, s, V, Bg, "std");
    bool svd_success = svd(U, s, V, Bg);

    if (svd_success == false) {
        pinv_success = false;
        inv_mat.zeros();
    } else {
        l = s.n_elem;
        arma::vec s1(l);
        tol = s(0)*1.490116e-08;
        for (i=0; i<l; i++){
            if (s(i) > tol) {
                s1(i) = s(i);
                r = r + 1;
            } else {
                s1(i) = 0.0;
            }
        }
        inv_mat = V.head_cols(r) * diagmat(1/s1.head(r)) * U.head_cols(r).t();
    }

    return Rcpp::List::create(pinv_success, inv_mat);
}



// [[Rcpp::export]]
Rcpp::List rrr_my(arma::mat X,
                  arma::mat Y,
                  int r){
    /*
     almost the same as RRR in R
     */
    bool rrr_success = true;
    Rcpp::List inv_X0;
    inv_X0 = pinv_my(X.t()*X);
    bool inv_X0_err;
    inv_X0_err = inv_X0[0];
    arma::mat Bl;
    arma::mat Al;
    arma::mat inv_X_mat;

    if (inv_X0_err == true){
        arma::mat inv_X_mat0 = inv_X0[1];
        inv_X_mat = inv_X_mat0;
    } else if (inv_X0_err == false) {
        int p = X.n_cols;
        Rcpp::List inv_X = pinv_my(X.t()*X + 0.1 * arma::eye(p, p));
        arma::mat inv_X_mat0 = inv_X[1];
        inv_X_mat = inv_X_mat0;
     }

    arma::mat fit_now = inv_X_mat * X.t() * Y;
    arma::mat proj_mat = X * fit_now;

    arma::mat U1;
    arma::vec s1;
    arma::mat V1;
    //svd(U1, s1, V1, proj_mat, "std");
    bool svd_success = true;
    svd_success = svd(U1, s1, V1, proj_mat);
    if (svd_success == true){
        Bl = fit_now * V1.cols(0, r-1);
        Al = V1.cols(0, r-1);
    } else {
        arma::mat SS = Y.t() * X * fit_now;
        SS = (SS + SS.t())/2.0;
        arma::vec eigval;
        arma::mat eigvec;

        rrr_success = eig_sym(eigval, eigvec, SS);

        arma::mat SS_r = eigvec.tail_cols(r);
        Bl = fit_now * SS_r;
        Al = SS_r;
    }
    return Rcpp::List::create(Bl, Al, rrr_success);
}


// [[Rcpp::export]]
Rcpp::List BLCD(arma::mat X,
                arma::mat Y,
                arma::mat Xl0,
                arma::mat Yl0,
                arma::mat Ag0,
                arma::mat Bg0,
                arma::mat Al0,
                arma::mat Bl0,
                int n,
                int r,
                int d,
                int p,
                int rx,
                int ry,
                int jx,
                int jy,
                int maxiter,
                double conv
                ){

    /* Notations:
     Ag -> U
     Bg -> V
     Bl0 -> A
     Bl0 -> B
     Xl, Xl0 -> X (I_jx \otimes V)
     Yl, Yl0 -> Y (I_jy \otimes U)
     */


    arma::vec obj(maxiter+1);
    obj.fill(0.0);
    arma::vec err_flag(maxiter+1);
    err_flag.zeros();
    double objnow;
    int iter = 0;
    int j, a1, a2, b1, b2;


    /* begin of Obj function */
    arma::mat Bl(jx*p, r);   // here Bl -> (I_jx \otimes V)*B
    arma::mat Al(jy*d, r);   // here Al -> (I_jy \otimes U)*A
    arma::mat C(jx*p, jy*d); // here C -> Y=XC
    arma::mat XC(n, jy*d);
    arma::mat E(n, jy*d);
    double sse;

    if ( Bg0.n_rows != Bg0.n_cols ){
        Bl = kron(arma::eye(jx, jx), Bg0) * Bl0;
    } else {
        Bl = Bl0;
    }
    if ( Ag0.n_rows != Ag0.n_cols ) {
        Al = kron(arma::eye(jy, jy), Ag0) * Al0;
    } else {
        Al = Al0;
    }

    C = Bl * Al.t();
    XC = X * C;
    E = Y - XC;
    sse = accu(E % E);
    /* end of Obj function */


    obj(0) = sse;
    objnow = sse + 10.0;


    arma::mat Ag1(d, ry);
    arma::mat Bg1(p, rx);
    arma::mat Al1(jy*ry, r);
    arma::mat Bl1(jx*rx, r);
    arma::mat Xg(n, jy*ry);
    arma::mat Yg0(d, ry);
    arma::vec yB(n*r);
    arma::mat XB(n * r, rx * p);
    arma::vec Bg1_vec(p*rx);
    //arma::mat fit_now(jx*rx, jy*ry);
    //arma::mat proj_mat(jy*ry, jy*ry);



    while ( (iter < maxiter) && (std::abs(obj(iter) - objnow) > conv) ) {

        /* update Ag (U) */
        if ( d == ry ) {
            Ag1 = arma::eye(ry, ry);
            Yl0 = Y;
        } else {
            Xg = Xl0 * Bl0 * Al0.t();
            Yg0.zeros();

            for (j = 0; j < jy; j++) {
                a1 = d * (j + 1 - 1) + 1 - 1;
                b1 = d * (j + 1) - 1;

                a2 = ry * (j + 1 - 1) + 1 - 1;
                b2 = ry * (j + 1) - 1;

                Yg0 = Yg0 + Y.cols(a1,b1).t() * Xg.cols(a2,b2);
            }
            arma::mat U;
            arma::vec s;
            arma::mat V;
            bool svd_Ag = svd_econ( U, s, V, Yg0, "both", "std" );
            if (svd_Ag == false) {
                err_flag(iter) = 1;
                break;
            } else {
                err_flag(iter) = 0;
            }
            Ag1 = U * V.t();
            Yl0 = Y * kron(arma::eye(jy, jy), Ag1);
        }


        /* update Bg (V) */
        if ( p == rx) {
            Bg1 = arma::eye(rx, rx);
            Xl0 = X;
        } else {
            yB = vectorise(Yl0 * Al0);
            XB.zeros();

            for (j = 0; j < jx; j++) {
                a1 = rx * (j + 1 - 1) + 1 - 1;
                b1 = rx * (j + 1) - 1;

                a2 = p * (j + 1 - 1) + 1 - 1;
                b2 = p * (j + 1) - 1;

                XB = XB + kron( Bl0.rows(a1, b1).t(), X.cols(a2, b2));
            }
            //Bg1_vec = pinv(XB.t()*XB + 0.0000001 * arma::eye(rx*p, rx*p), 1.490116e-08, "std") * XB.t() * yB;
            Rcpp::List inv_XB = pinv_my(XB.t()*XB + 0.0000001 * arma::eye(rx*p, rx*p));
            bool inv_XB_err = inv_XB[0];
            arma::mat inv_XB_mat = inv_XB[1];

            if (inv_XB_err == false) {
                err_flag(iter) = 2;
                break;
            }

            Bg1_vec = inv_XB_mat * XB.t() * yB;
            Bg1 = reshape(Bg1_vec, p, rx);

            arma::mat Q;
            arma::mat R;

            bool svd_Bg = arma::qr_econ(Q, R, Bg1);
            if (svd_Bg == false) {
                err_flag(iter) = 3;
                break;
            }
            Bg1 = Q;

            Xl0 = X * kron(arma::eye(jx, jx), Bg1);
        }

        /* update Al (A) and Bl (B) */

        /* RRR start */
        Rcpp::List inv_X0 = pinv_my(Xl0.t()*Xl0);
        bool inv_X0_err = inv_X0[0];

        if (inv_X0_err == true){
            arma::mat inv_X_mat = inv_X0[1];

            arma::mat fit_now = inv_X_mat * Xl0.t() * Yl0;
            arma::mat proj_mat = Xl0 * fit_now;

            arma::mat U2;
            arma::vec s2;
            arma::mat V2;

            bool svd_RRR = svd(U2, s2, V2, proj_mat, "std");
            if (svd_RRR == false) {
                err_flag(iter) = 4;
                break;
            }

            Bl1 = fit_now * V2.cols(0, r-1);
            Al1 = V2.cols(0, r-1);

        } else if (inv_X0_err == false) {
            Rcpp::List inv_X = pinv_my(Xl0.t()*Xl0 + 0.1 * arma::eye(jx *rx, jx*rx));
            arma::mat inv_X_mat = inv_X[1];

            arma::mat fit_now = inv_X_mat * Xl0.t() * Yl0;
            arma::mat proj_mat = Xl0 * fit_now;

            arma::mat U2;
            arma::vec s2;
            arma::mat V2;

            bool svd_RRR = svd(U2, s2, V2, proj_mat, "std");
            if (svd_RRR == false) {
                err_flag(iter) = 4;
                break;
            }

            Bl1 = fit_now * V2.cols(0, r-1);
            Al1 = V2.cols(0, r-1);
        }
        /* RRR end */



        iter++;

        /* start of Obj function */
        if ( Bg1.n_rows != Bg1.n_cols ){
            Bl = kron(arma::eye(jx, jx), Bg1) * Bl1;
        } else {
            Bl = Bl1;
        }
        if ( Ag1.n_rows != Ag1.n_cols ) {
            Al = kron(arma::eye(jy, jy), Ag1) * Al1;
        } else {
            Al = Al1;
        }

        C = Bl * Al.t();
        XC = X * C;
        E = Y - XC;
        sse = accu(E % E);
        /* end of Obj function */

        obj(iter) = sse;
        objnow = obj(iter - 1);

        Al0 = Al1;
        Bl0 = Bl1;
        Ag0 = Ag1;
        Bg0 = Bg1;
    }

    //formulate output
    arma::vec output_obj = obj.subvec(0, iter);
    arma::vec output_err = err_flag;
    int output_iter = iter + 1;

    return Rcpp::List::create(Rcpp::Named("sse") = sse,
                              Rcpp::Named("C") = C,
                              Rcpp::Named("Ag1") = Ag1,
                              Rcpp::Named("Bg1") = Bg1,
                              Rcpp::Named("Al1") = Al0,
                              Rcpp::Named("Bl1") = Bl0,
                              Rcpp::Named("obj") = output_obj,
                              Rcpp::Named("err") = output_err,
                              Rcpp::Named("iter") = output_iter
                              );}




// [[Rcpp::export]]
Rcpp::List nrrr_init_my(
                        arma::mat Y,
                        arma::mat X,
                        int r,
                        int rx,
                        int ry,
                        int jx,
                        int jy,
                        int p,
                        int d,
                        int n
                        ){
    Rcpp::List Cr_list = rrr_my(X, Y, r);
    arma::mat Cr1 = Cr_list[0];
    arma::mat Cr2 = Cr_list[1];
    arma::mat Crr = Cr1 * Cr2.t();

    arma::mat Bghat;
    arma::mat Aghat;
    arma::mat Cbg(p, d*jy*jx);
    arma::mat Cag(d, p*jy*jx);
    int j, a1, b1, a2, b2;

    /* get V */
    if ( p == rx ){
        Bghat = arma::eye(rx, rx);
    } else {
        for (j=0; j<jx; j++){
            a1 = d * jy * j;
            b1 = d * jy * (j + 1) - 1;

            a2 = p * j;
            b2 = p * (j + 1) - 1;

            Cbg.cols(a1, b1) = Crr.rows(a2, b2);
        }
        arma::mat U;
        arma::vec s;
        arma::mat V;

        //svd(U, s, V, Cbg, "std");
        svd(U, s, V, Cbg);
        Bghat = U.cols(0, rx-1);
    }

    /* get U */
    if ( d == ry ){
        Aghat = arma::eye(ry, ry);
    } else {
        for (j=0; j<jy; j++){
            a1 = p * jx * j;
            b1 = p * jx * (j + 1) - 1;

            a2 = d * j;
            b2 = d * (j + 1) - 1;

            Cag.cols(a1, b1) = Crr.cols(a2, b2).t();
        }
        arma::mat U1;
        arma::vec s1;
        arma::mat V1;

        //svd(U1, s1, V1, Cag, "std");
        svd(U1, s1, V1, Cag);
        Aghat = U1.cols(0, ry-1);
    }

    return Rcpp::List::create(Aghat, Bghat, Crr);
}



// [[Rcpp::export]]
Rcpp::List nrrr_est_my(
                arma::mat Y,
                arma::mat X,
                int rini,
                int r,
                int rx,
                int ry,
                int jx,
                int jy,
                int p,
                int d,
                int n,
                int maxiter,
                double conv,
                int method,
                double lambda
                ){

    /* method:  1-RRR, 2-RRS */

    if ( ( method == 2 )&&( lambda > 0.0 ) ){
        Y = join_cols(Y, arma::zeros(p*jx, d*jy));
        X = join_cols(X, sqrt(lambda)*arma::eye(p*jx, p*jx));
        n = n + p*jx;
    }

    arma::mat Ag0;
    arma::mat Bg0;

    //if ( (Ag0 == NULL) || (Bg0 == NULL) ){
        Rcpp::List ini_res = nrrr_init_my(Y, X, rini, rx, ry, jx, jy, p, d, n);
        arma::mat ini_res_Ag0 = ini_res[0];
        arma::mat ini_res_Bg0 = ini_res[1];

    Ag0 = ini_res_Ag0;
    Bg0 = ini_res_Bg0;
    //}

    arma::mat Yl0;
    arma::mat Xl0;

    if (d == ry){
        Yl0 = Y;
    } else {
        Yl0 = Y * kron(arma::eye(jy, jy), Ag0);
    }

    if (p == rx){
        Xl0 = X;
    } else {
        Xl0 = X * kron(arma::eye(jx, jx), Bg0);
    }

    /* given U and V, compute Al, Bl */
    arma::mat Bl0;
    arma::mat Al0;

    Rcpp::List fitRR = rrr_my(Xl0, Yl0, r);
    arma::mat fitRR_Bl0 = fitRR[0];
    arma::mat fitRR_Al0 = fitRR[1];
    Bl0 = fitRR_Bl0;
    Al0 = fitRR_Al0;


    /*---------------------- BLCD start ----------------------*/
    arma::vec obj(maxiter+1);
    obj.fill(0.0);
    arma::vec err_flag(maxiter+1);
    err_flag.zeros();
    int any_err = 0;
    double objnow;
    int iter = 0;
    int j, a1, a2, b1, b2;


    /*-------- begin of Obj function -------*/
    arma::mat Bl(jx*p, r);   // here Bl -> (I_jx \otimes V)*B
    arma::mat Al(jy*d, r);   // here Al -> (I_jy \otimes U)*A
    arma::mat C(jx*p, jy*d); // here C -> Y=XC
    arma::mat XC(n, jy*d);
    arma::mat E(n, jy*d);
    double sse;

    if ( Bg0.n_rows != Bg0.n_cols ){
        Bl = kron(arma::eye(jx, jx), Bg0) * Bl0;
    } else {
        Bl = Bl0;
    }
    if ( Ag0.n_rows != Ag0.n_cols ) {
        Al = kron(arma::eye(jy, jy), Ag0) * Al0;
    } else {
        Al = Al0;
    }

    C = Bl * Al.t();
    XC = X * C;
    E = Y - XC;
    sse = accu(E % E);
    /* ------- end of Obj function ---------*/


    obj(0) = sse;
    objnow = sse + 10.0;


    arma::mat Ag1(d, ry);
    arma::mat Bg1(p, rx);
    arma::mat Al1(jy*ry, r);
    arma::mat Bl1(jx*rx, r);
    arma::mat Xg(n, jy*ry);
    arma::mat Yg0(d, ry);
    arma::vec yB(n*r);
    arma::mat XB(n * r, rx * p);
    arma::vec Bg1_vec(p*rx);


    while ( (iter < maxiter) && (std::abs(obj(iter) - objnow) > conv) ) {

        /* update Ag (U) */
        if ( d == ry ) {
            Ag1 = arma::eye(ry, ry);
            Yl0 = Y;
        } else {
            Xg = Xl0 * Bl0 * Al0.t();
            Yg0.zeros();

            for (j = 0; j < jy; j++) {
                a1 = d * j;
                b1 = d * (j + 1) - 1;

                a2 = ry * j;
                b2 = ry * (j + 1) - 1;

                Yg0 = Yg0 + Y.cols(a1,b1).t() * Xg.cols(a2,b2);
            }
            arma::mat U;
            arma::vec s;
            arma::mat V;
            //bool svd_Ag = svd_econ( U, s, V, Yg0, "both", "std" );
            bool svd_Ag = svd_econ( U, s, V, Yg0, "both");
            if (svd_Ag == false) {
                err_flag(iter) = 1;
                break;
            } else {
                err_flag(iter) = 0;
            }
            Ag1 = U * V.t();
            Yl0 = Y * kron(arma::eye(jy, jy), Ag1);
        }


        /* update Bg (V) */
        if ( p == rx) {
            Bg1 = arma::eye(rx, rx);
            Xl0 = X;
        } else {
            yB = vectorise(Yl0 * Al0);
            XB.zeros();

            for (j = 0; j < jx; j++) {
                a1 = rx * j;
                b1 = rx * (j + 1) - 1;

                a2 = p * j;
                b2 = p * (j + 1) - 1;

                XB = XB + kron( Bl0.rows(a1, b1).t(), X.cols(a2, b2));
            }
            Rcpp::List inv_XB = pinv_my(XB.t()*XB + 0.0000001 * arma::eye(rx*p, rx*p));
            bool inv_XB_err = inv_XB[0];
            arma::mat inv_XB_mat = inv_XB[1];

            if (inv_XB_err == false) {
                err_flag(iter) = 2;
                break;
            }

            Bg1_vec = inv_XB_mat * XB.t() * yB;
            Bg1 = reshape(Bg1_vec, p, rx);

            arma::mat Q;
            arma::mat R;

            bool svd_Bg = arma::qr_econ(Q, R, Bg1);
            if (svd_Bg == false) {
                err_flag(iter) = 3;
                break;
            }
            Bg1 = Q;

            Xl0 = X * kron(arma::eye(jx, jx), Bg1);
        }

        /* update Al (A) and Bl (B) */

        Rcpp::List fitRR = rrr_my(Xl0, Yl0, r);
        arma::mat fitRR_Bl1 = fitRR[0];
        arma::mat fitRR_Al1 = fitRR[1];
        bool fitRR_success = fitRR[2];
        if (fitRR_success == false) {
            err_flag(iter) = 4;
            break;
        }
        Bl1 = fitRR_Bl1;
        Al1 = fitRR_Al1;

        /* --------- RRR start ----------
        Rcpp::List inv_X0 = pinv_my(Xl0.t()*Xl0);
        bool inv_X0_err = inv_X0[0];

        if (inv_X0_err == true){
            arma::mat inv_X_mat = inv_X0[1];

            arma::mat fit_now = inv_X_mat * Xl0.t() * Yl0;
            arma::mat proj_mat = Xl0 * fit_now;

            arma::mat U2;
            arma::vec s2;
            arma::mat V2;

            //bool svd_RRR = svd(U2, s2, V2, proj_mat, "std");
            bool svd_RRR = svd(U2, s2, V2, proj_mat);
            if (svd_RRR == false) {
                err_flag(iter) = 4;
                break;
            }

            Bl1 = fit_now * V2.cols(0, r-1);
            Al1 = V2.cols(0, r-1);

        } else if (inv_X0_err == false) {
            Rcpp::List inv_X = pinv_my(Xl0.t()*Xl0 + 0.1 * arma::eye(jx *rx, jx*rx));
            arma::mat inv_X_mat = inv_X[1];

            arma::mat fit_now = inv_X_mat * Xl0.t() * Yl0;
            arma::mat proj_mat = Xl0 * fit_now;

            arma::mat U2;
            arma::vec s2;
            arma::mat V2;

            //bool svd_RRR = svd(U2, s2, V2, proj_mat, "std");
            bool svd_RRR = svd(U2, s2, V2, proj_mat);
            if (svd_RRR == false) {
                err_flag(iter) = 4;
                break;
            }

            Bl1 = fit_now * V2.cols(0, r-1);
            Al1 = V2.cols(0, r-1);
        }
        -------- RRR end ---------*/



        iter++;

        /*-------- start of Obj function --------*/
        if ( Bg1.n_rows != Bg1.n_cols ){
            Bl = kron(arma::eye(jx, jx), Bg1) * Bl1;
        } else {
            Bl = Bl1;
        }
        if ( Ag1.n_rows != Ag1.n_cols ) {
            Al = kron(arma::eye(jy, jy), Ag1) * Al1;
        } else {
            Al = Al1;
        }

        C = Bl * Al.t();
        XC = X * C;
        E = Y - XC;
        sse = accu(E % E);
        /*-------- end of Obj function --------*/

        obj(iter) = sse;
        objnow = obj(iter - 1);

        Al0 = Al1;
        Bl0 = Bl1;
        Ag0 = Ag1;
        Bg0 = Bg1;
    }
    if ( (iter == maxiter) && (std::abs(obj(iter) - objnow) > conv) ) {
        err_flag(iter) = 5;
    }
    /*---------------------- BLCD end ----------------------*/

    int err_which = sum(err_flag);
    if ( err_which != 0 ) {
        any_err = 1;
    }

    double xr = arma::rank(X, 0.01);
    double df0;
    if (xr > jx*rx) {
        df0 = rx * (xr/jx - rx);
    } else {
        df0 = 0.0;
    }
    double df = df0 + ry * (d - ry) + (jy * ry + jx * rx - r) * r;

    if (sse < 0.1) {
        sse = 0.0;
    }
    if (method == 2) {
        sse = sse - lambda * accu(C % C);
    }

    arma::vec ic(4);
    ic(0) = log(sse) + log(d * jy * n) * df / (d * jy * n) ;
    ic(1) = log(sse) + 2 * log(d * jy * p * jx) * df / (d * jy * n);
    ic(2) = log(sse) + 2 * df / (d * jy * n);
    ic(3) = sse / (d * jy * n * (1 - df / (d * jy * n)) * (1 - df / (d * jy * n)));


    //formulate output
    arma::vec output_obj = obj.subvec(0, iter);
    int output_iter = iter + 1;


    return Rcpp::List::create(Rcpp::Named("Ag") = Ag1,
                              Rcpp::Named("Bg") = Bg1,
                              Rcpp::Named("Al") = Al1,
                              Rcpp::Named("Bl") = Bl1,
                              Rcpp::Named("C") = C,
                              Rcpp::Named("df") = df,
                              Rcpp::Named("err_flag") = any_err,
                              Rcpp::Named("err_which") = err_which,
                              Rcpp::Named("sse") = sse,
                              Rcpp::Named("ic") = ic,
                              Rcpp::Named("obj") = output_obj,
                              Rcpp::Named("iter") = output_iter
                              );
}



// [[Rcpp::export]]
Rcpp::List nrrr_select_my(arma::mat Y,
                          arma::mat X,
                          int xr,
                          int rfit,
                          int jx,
                          int jy,
                          int p,
                          int d,
                          int n,
                          int ic,
                          int maxiter,
                          double conv,
                          int method,
                          double lambda,
                          int dimred1,
                          int dimred2,
                          int dimred3
                          ){

    /* ic: 0-BIC, 1-BICP, 2-AIC, 3-GCV */
    /* method: 1-RRR, 2-RRS            */

    int rest = rfit;
    int i, j, l;
    int rxfit, ryfit;
    int rxest, ryest;
    int r_length = std::min(jx*p, jy*d);
    arma::uword min_ind;
    arma::vec rxErrseq(p);
    arma::vec ryErrseq(d);
    arma::vec rErrseq(r_length);

    /* select rx */
    if ( dimred2 == 1 ){
        arma::vec rxfitseq(p);
        arma::vec icseq_x(p);
        ryfit = d;

        for (i=0; i<p; i++){
            rxfitseq(i) = i + 1;
            rxfit = i + 1;

            Rcpp::List fit1 = nrrr_est_my(Y,X,rfit,rfit,rxfit,ryfit,jx,jy,p,
                                          d,n,maxiter,conv,
                                          method,lambda);
            arma::vec fit1_ic = fit1[9];
            icseq_x(i) = fit1_ic(ic);
            rxErrseq(i) = fit1[7];
        }
        min_ind = index_min(icseq_x);
        rxest = rxfitseq(min_ind);
    } else {
        rxest = p;
    }
    rxfit = rxest;


    /* select ry */
    if ( dimred3 == 1 ){
        arma::vec ryfitseq(d);
        arma::vec icseq_y(d);

        for (i=0; i<d; i++){
            ryfitseq(i) = i + 1;
            ryfit = i + 1;

            Rcpp::List fit1 = nrrr_est_my(Y,X,rfit,rfit,rxfit,ryfit,jx,jy,p,
                                          d,n,maxiter,conv,
                                          method,lambda);
            arma::vec fit1_ic = fit1[9];
            icseq_y(i) = fit1_ic(ic);
            ryErrseq(i) = fit1[7];
        }
        min_ind = index_min(icseq_y);
        ryest = ryfitseq(min_ind);
    } else {
        ryest = d;
    }
    ryfit = ryest;


    /* refine r */
    if ( dimred1 == 1 ){

        arma::vec r_range(4);
        r_range(0) = n;
        r_range(1) = p*jx;
        r_range(2) = d*jy;
        r_range(3) = rest + 5;

        int a = std::max(1, rest - 5);
        int b = arma::min(r_range);

        r_length = b - a + 1;
        arma::vec rfitseq(r_length);

        for (j=a; j< b + 1; j++){
            rfitseq(j-a) = j;
        }

        arma::vec icseq_r(r_length);

        for (l=0; l< r_length; l++){
            rfit = rfitseq(l);

            Rcpp::List fit1 = nrrr_est_my(Y,X,rfit,rfit,rxfit,ryfit,jx,jy,p,
                                          d,n,maxiter,conv,
                                          method,lambda);
            arma::vec fit1_ic = fit1[9];
            icseq_r(l) = fit1_ic(ic);
            rErrseq(l) = fit1[7];
        }
        min_ind = index_min(icseq_r);
        rest = rfitseq(min_ind);
    }


    Rcpp::List fit_final = nrrr_est_my(Y,X,rest,rest,rxfit,ryfit,jx,jy,p,
                                  d,n,maxiter,conv,
                                  method,lambda);

    arma::vec rErrseq_output = rErrseq.subvec(0, r_length-1);


    return Rcpp::List::create(Rcpp::Named("Ag") = fit_final[0],
                              Rcpp::Named("Bg") = fit_final[1],
                              Rcpp::Named("Al") = fit_final[2],
                              Rcpp::Named("Bl") = fit_final[3],
                              Rcpp::Named("C") = fit_final[4],
                              Rcpp::Named("df") = fit_final[5],
                              Rcpp::Named("sse") = fit_final[8],
                              Rcpp::Named("ic") = fit_final[9],
                              Rcpp::Named("obj") = fit_final[10],
                              Rcpp::Named("rank") = rest,
                              Rcpp::Named("rx") = rxest,
                              Rcpp::Named("ry") = ryest,
                              Rcpp::Named("rxErrseq") = rxErrseq,
                              Rcpp::Named("ryErrseq") = ryErrseq,
                              Rcpp::Named("rErrseq") = rErrseq_output
                              );
}


// [[Rcpp::export]]
arma::mat del_rows(arma::mat X, arma::uvec e){
    X.shed_rows(e);
    return X;
}


// [[Rcpp::export]]
Rcpp::List nrrr_cv_my(arma::mat Y,
                      arma::mat X,
                      arma::uvec norder,
                      int nfold,
                      int xr,
                      int rfit,
                      int xrankfix,
                      int yrankfix,
                      int jx,
                      int jy,
                      int p,
                      int d,
                      int n,
                      int maxiter,
                      int method,
                      int dimred1,
                      int dimred2,
                      int dimred3,
                      double conv,
                      double lambda
                      ){

    /* method: 1-RRR, 2-RRS                  */
    /* xrankfix: 0-null, others-specified rx */

    int rest;
    rest = rfit;

    int ndel;
    int f, i, nf;
    arma::uvec iddel;
    arma::mat Xf;
    arma::mat Xfdel;
    arma::mat Yf;
    arma::mat Yfdel;
    int rxfit, ryfit;
    int rxest, ryest;
    arma::mat rxErrmat(p,nfold);
    rxErrmat.zeros();
    arma::mat ryErrmat(d,nfold);
    ryErrmat.zeros();

    ndel = round(n/nfold);

    /* select rx */
    arma::mat rx_path(p, nfold);
    rx_path.zeros();
    if (p == 1) {
        rxest = p;
    } else {
        if (dimred2 == 1){
            arma::vec rxfitseq(p);
            for (f=0; f<nfold; f++){
                if (f != nfold - 1){
                    iddel = norder.subvec(ndel * f, ndel * (f + 1) - 1);
                } else {
                    iddel = norder.subvec(ndel * f, n-1);
                }
                ndel = iddel.n_elem;
                nf = n - ndel;

                Xf = del_rows(X, iddel);
                Xfdel = X.rows(iddel);
                Yf = del_rows(Y, iddel);
                Yfdel = Y.rows(iddel);

                ryfit = d;
                for (i=0; i<p; i++){
                    rxfit = i + 1;
                    rxfitseq(i) = i + 1;

                    Rcpp::List fit1 = nrrr_est_my(Yf,Xf,rfit,rfit,rxfit,ryfit,jx,jy,p,
                                                  d,nf,maxiter,conv,
                                                  method,lambda);
                    arma::mat C_est = fit1[4];
                    rx_path(i, f) = accu( arma::square(Yfdel - Xfdel * C_est) );
                    rxErrmat(i, f) = fit1[7];
                }
            }
            arma::vec crerr = arma::sum(rx_path, 1);
            rxest = rxfitseq(index_min(crerr));
        } else {
            if (xrankfix == 0) {
                rxest = p;
            } else {
                rxest = xrankfix;
            }
        }
    }
    rxfit = rxest;


    /* select ry */
    arma::mat ry_path(d, nfold);
    ry_path.zeros();
    if (d == 1) {
        ryest = d;
    } else {
        if (dimred3 == 1){

            arma::vec ryfitseq(d);
            for (f=0; f<nfold; f++){
                if (f != nfold - 1){
                    iddel = norder.subvec(ndel * f, ndel * (f + 1) - 1);
                } else {
                    iddel = norder.subvec(ndel * f, n-1);
                }
                ndel = iddel.n_elem;
                nf = n - ndel;

                Xf = del_rows(X, iddel);
                Xfdel = X.rows(iddel);
                Yf = del_rows(Y, iddel);
                Yfdel = Y.rows(iddel);

                for (i=0; i<d; i++){
                    ryfit = i + 1;
                    ryfitseq(i) = i + 1;

                    Rcpp::List fit1 = nrrr_est_my(Yf,Xf,rfit,rfit,rxfit,ryfit,jx,jy,p,
                                                  d,nf,maxiter,conv,
                                                  method,lambda);
                    arma::mat C_est = fit1[4];
                    ry_path(i, f) = accu( arma::square(Yfdel - Xfdel * C_est) );
                    ryErrmat(i, f) = fit1[7];
                }
            }
            arma::vec crerr = arma::sum(ry_path, 1);
            ryest = ryfitseq(index_min(crerr));
        } else {
            if (yrankfix == 0) {
                ryest = d;
            } else {
                ryest = yrankfix;
            }
        }
    }
    ryfit = ryest;

    /* refine r */
    arma::vec r_range(4);
    r_range(0) = n;
    r_range(1) = p*jx;
    r_range(2) = d*jy;
    r_range(3) = rest + 5;

    int a = std::max(1, rest - 5);
    int b = arma::min(r_range);

    int r_length = b - a + 1;
    arma::vec rfitseq(r_length);

    arma::mat r_path(r_length, nfold);
    r_path.zeros();
    arma::mat rErrmat(r_length, nfold);
    rErrmat.zeros();
    if ( dimred1 == 1 ){

        for (f=0; f<nfold; f++){
            if (f != nfold - 1){
                iddel = norder.subvec(ndel * f, ndel * (f + 1) - 1);
            } else {
                iddel = norder.subvec(ndel * f, n-1);
            }
            ndel = iddel.n_elem;
            nf = n - ndel;

            Xf = del_rows(X, iddel);
            Xfdel = X.rows(iddel);
            Yf = del_rows(Y, iddel);
            Yfdel = Y.rows(iddel);

            for (i=a; i<b+1; i++){
                rfit = i;
                rfitseq(i-a) = i;

                Rcpp::List fit1 = nrrr_est_my(Yf,Xf,rfit,rfit,rxfit,ryfit,jx,jy,p,
                                              d,nf,maxiter,conv,
                                              method,lambda);
                arma::mat C_est = fit1[4];
                r_path(i-a, f) = accu( arma::square(Yfdel - Xfdel * C_est) );
                rErrmat(i-a, f) = fit1[7];
            }
        }
        arma::vec crerr = arma::sum(r_path, 1);
        rest = rfitseq(index_min(crerr));
    }

    Rcpp::List fit_final = nrrr_est_my(Y,X,rest,rest,rxfit,ryfit,jx,jy,p,
                                      d,n,maxiter,conv,
                                      method,lambda);

    return Rcpp::List::create(Rcpp::Named("Ag") = fit_final[0],
                              Rcpp::Named("Bg") = fit_final[1],
                              Rcpp::Named("Al") = fit_final[2],
                              Rcpp::Named("Bl") = fit_final[3],
                              Rcpp::Named("C") = fit_final[4],
                              Rcpp::Named("df") = fit_final[5],
                              Rcpp::Named("sse") = fit_final[8],
                              Rcpp::Named("ic") = fit_final[9],
                              Rcpp::Named("obj") = fit_final[10],
                              Rcpp::Named("rx_path") = rx_path,
                              Rcpp::Named("ry_path") = ry_path,
                              Rcpp::Named("r_path") = r_path,
                              Rcpp::Named("rank") = rest,
                              Rcpp::Named("rx") = rxest,
                              Rcpp::Named("ry") = ryest,
                              Rcpp::Named("rxErrmat") = rxErrmat,
                              Rcpp::Named("ryErrmat") = ryErrmat,
                              Rcpp::Named("rErrmat") = rErrmat
                              );
}
