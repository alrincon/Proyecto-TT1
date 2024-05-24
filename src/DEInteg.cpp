#include "../include/DEInteg.h"

Matrix extractColt(Matrix* input, int row){
    Matrix res(1, input->getColumnas());
    for(int i = 1; i <= input->getColumnas(); i++){
        res(1, i) = (*input)(row, 1);
    }

    return res;
}

struct DE_STATE {
    static const int DE_INIT = 1;       // Restart integration
    static const int DE_DONE = 2;       // Successful step
    static const int DE_BADACC = 3;     // Accuracy requirement could not be achieved
    static const int DE_NUMSTEPS = 4;   // Permitted number of steps exceeded
    static const int DE_STIFF = 5;      // Stiff problem suspected
    static const int DE_INVPARAM = 6;   // Invalid input parameters
};

Matrix DEInteg(Matrix (*func)(double, Matrix*), double t, double tout, double relerr, double abserr, int n_eqn, Matrix *yInput) {
    Matrix y(yInput->getFilas(), yInput->getColumnas());
    y.redefine(yInput);

    double twou = 2 * 2.2204e-16;
    double fouru = 4 * 2.2204e-16;

    int limitn = 100;
    int State_ = DE_STATE::DE_INIT;
    bool PermitTOUT = true; // Allow integration past tout by default
    bool OldPermit;
    double told = 0;

    Matrix two(1,14);
    for (int i = 1; i <= 14; i++) {
        two(1,i) = pow(2.0, i-1);
    }

    Matrix gstr(1,14);
    gstr(1, 1) = 1.0;
    gstr(1,2) = 0.5;
    gstr(1,3) = 0.0833;
    gstr(1,4) = 0.0417;
    gstr(1,5) = 0.0264;
    gstr(1,6) = 0.0188;
    gstr(1,7) = 0.0143;
    gstr(1,8) = 0.0114;
    gstr(1,9) = 0.00936;
    gstr(1,10) = 0.00789;
    gstr(1,11) = 0.00679;
    gstr(1,12) = 0.00592;
    gstr(1,13) = 0.00524;
    gstr(1,14) = 0.00468;

    //Hay que darle la vuelta a todos
    Matrix yy(1,n_eqn);    // Allocate vectors with proper dimension
    Matrix wt(1,n_eqn);
    Matrix p(1,n_eqn);
    Matrix yp(1,n_eqn);
    Matrix phi(17,n_eqn);
    Matrix g( 1,14);
    Matrix sig( 1,14);
    Matrix rho( 1,14);
    Matrix w( 1,13);
    Matrix alpha(1,13);
    Matrix beta(1,13);
    Matrix v(1,13);
    Matrix psi_(1,13);

    if (t == tout) return y;

    double epsilon = fmax(relerr, abserr);
    if (relerr < 0.0 || abserr < 0.0 || epsilon <= 0.0 || State_ > DE_STATE::DE_INVPARAM || (State_ != DE_STATE::DE_INIT && t != told)) {
        State_ = DE_STATE::DE_INVPARAM;
        cout << "SALIDA DE_INVPARAM" << endl;
        return y;
    }

    double del = tout - t;
    double absdel = fabs(del);

    double tend = t + 100.0 * del;
    if (!PermitTOUT) tend = tout;

    int nostep = 0;
    int nnostep = 0;
    int kle4 = 0;
    bool stiff = false;
    double releps = relerr / epsilon;
    double abseps = abserr / epsilon;

    bool start = true;
    double x=0;
    double delsgn=0;
    double h=0;


    if  ( (State_==DE_STATE::DE_INIT) || (delsgn*del<=0.0) ){
        // On start and restart also set the work variables x and yy(*),
        // store the direction of integration and initialize the step size
        start  = true;
        x      = t;
        yy     = (y);
        delsgn = sign(1.0, del);
        h      = sign( fmax(fouru*abs(x), fabs(tout-x)), tout-x );
    }

    Matrix yout(1,n_eqn);
    Matrix ypout(1,n_eqn);
    Matrix typ(1,n_eqn);
    while (true) {
        double hi;
        int ki;
        int kold = 0;
        nnostep++;

        // If already past output point, interpolate solution and return

        if (fabs(x - t) >= absdel) {
            g( 1,2) = 1.0;
            rho( 1,2) = 1.0;
            hi = tout - x;
            ki = kold + 1;

            // Initialize w[*] for computing g[*]
            for (int i = 1; i <= ki; i++) {
                double temp1 = i;
                w(1,i + 1) = 1.0 / temp1;
            }

            // Compute g[*]
            double term = 0.0;

            for (int j = 2; j <= ki; j++) {
                double psijm1 = psi_(1,j);
                double gamma = (hi + term) / psijm1;
                double eta = hi / psijm1;

                for (int i = 1; i <= ki + 1 - j; i++) {
                    w(1,i + 1) = gamma * w(1,i + 1) - eta * w(1,i + 2);
                }

                g(1,j + 1) = w(1, 2);
                rho(1,j + 1) = gamma * rho(1,j);
                term = psijm1;
            }


            // Interpolate for the solution yout and for
            // the derivative of the solution ypout
            for (int j = 1; j <= ki; j++) {
                int i = ki + 1 - j;
                yout = yout + extractColt(&phi, i + 1) * g(1,i + 1);
                ypout = ypout + extractColt(&phi, i + 1) * rho(1,i + 1);
            }

            yout = y + yout * hi;
            y = yout;

            State_ = DE_STATE::DE_DONE; // Set return code
            t = tout;             // Set independent variable
            told = t;                // Store independent variable
            OldPermit = PermitTOUT;
            cout << "SALIDA DE_DONE1" << endl;
            return y;                       // Normal exit
        }

        // If cannot go past output point and sufficiently close,
        // extrapolate and return
        if (!PermitTOUT && (fabs(tout - x) < fouru * fabs(x))) {
            h = tout - x;
            yp = func(x, &yy);          // Compute derivative yp(x)
            typ(1,n_eqn);
            typ = yy + (yp) * h;
            y = typ;                // Extrapolate vector from x to tout
            State_ = DE_STATE::DE_DONE; // Set return code
            t = tout;             // Set independent variable
            told = t;                // Store independent variable
            OldPermit = PermitTOUT;
            cout << "SALIDA DE_DONE2" << endl;
            return y;                       // Normal exit
        }

        // Limit step size, set weight vector and take a step
        if(nnostep > limitn){return y;}
        h = sign(fmin(fabs(h), fabs(tend - x)), h);
        for (int l = 1; l <= n_eqn; l++) {
            wt(1,l) = releps * fabs(yy(1,l)) + abseps;
        }

        //   Step
        //
        // Begin block 0
        //
        // Check if step size or error tolerance is too small for machine
        // precision.  If first step, initialize phi array and estimate a
        // starting step size. If step size is too small, determine an
        // acceptable one.
        //

        bool crash;

        if (fabs(h) < fouru * fabs(x)) {
            h = sign(fouru * fabs(x), h);
            crash = true;
            cout << "SALIDA crash1" << endl;
            return y;           // Exit
        }

        double p5eps = 0.5 * epsilon;
        crash = false;
        g(1,2) = 1.0;
        g(1,3) = 0.5;
        sig(1,2) = 1.0;

        double ifail = 0;

        // If error tolerance is too small, increase it to an
        // acceptable value.

        double round = 0.0;

        for (int l = 1; l <= n_eqn; l++) {
            round = round + ((y)(1,l) * (y)(1,l)) / (wt(1,l) * wt(1,l));
        }

        round = twou * sqrt(round);

        if (p5eps < round) {
            epsilon = 2.0 * round * (1.0 + fouru);
            crash = true;
            cout << "SALIDA crash2" << endl;
            return y;
        }


        double absh;
        int k = 0;
        double hold;
        double hnew;
        bool phase1;
        bool nornd;

        if (start) {
            // Initialize. Compute appropriate step size for first step.


            yp = func(x, &y);
            double sum = 0.0;

            for (int l = 1; l <= n_eqn; l++) {
                phi(2,l) = yp(1,l);
                phi(3,l) = 0.0;
                sum = sum + (yp(1,l) * yp(1,l)) / (wt(1,l) * wt(1,l));
            }

            sum = sqrt(sum);
            absh = fabs(h);

            if (epsilon < 16.0 * sum * h * h) {
                absh = 0.25 * sqrt(epsilon / sum);
            }

            double ttt = sign(fmax(absh, fouru * fabs(x)), h);
            h = ttt;
            hold = 0.0;
            hnew = 0.0;
            k = 1;
            kold = 0;
            start = false;
            phase1 = true;
            nornd = true;

            if (p5eps <= 100.0 * round) {
                nornd = false;
                for (int l = 1; l <= n_eqn; l++) {
                    phi(16,l) = 0.0;
                }
            }
        }

        //
        // End block 0
        //

        //
        // Repeat blocks 1, 2 (and 3) until step is successful
        //

        int kp1;
        int kp2;
        int km1;
        int km2;
        int knew;

        double erkm2;
        double erkm1;
        double erk;

        int ns;
        int nsp1;
        int realns;

        while (true) {
            //
            // Begin block 1
            //
            // Compute coefficients of formulas for this step. Avoid computing
            // those quantities not changed when step size is not changed.
            //

            kp1 = k + 1;
            kp2 = k + 2;
            km1 = k - 1;
            km2 = k - 2;

            // ns is the number of steps taken with size h, including the
            // current one. When k<ns, no coefficients change.

            double temp1;

            if (h != hold) {
                ns = 0;
            }

            if (ns <= kold) {
                ns = ns + 1;
            }

            nsp1 = ns + 1;

            if (k >= ns) {
                // Compute those components of alpha[*],beta[*],psi[*],sig[*]
                // which are changed
                beta(1,ns + 1) = 1.0;
                realns = ns;
                alpha(1,ns + 1) = 1.0 / realns;
                temp1 = h * realns;
                sig(1,nsp1 + 1) = 1.0;

                if (k >= nsp1) {
                    for (int i = nsp1; i <= k; i++) {
                        int im1 = i - 1;
                        double temp2 = psi_(1,im1 + 1);
                        psi_(1,im1 + 1) = temp1;
                        beta(1,i + 1) = beta(1,im1 + 1) * psi_(1,im1 + 1) / temp2;
                        temp1 = temp2 + h;
                        alpha(1,i + 1) = h / temp1;
                        int reali = i;
                        sig(1,i + 2) = reali * alpha(1,i + 1) * sig(1,i + 1);
                    }
                }

                psi_(1,k + 1) = temp1;

                // Compute coefficients g[*]; initialize v[*] and set w[*].
                if (ns > 1) {
                    // If order was raised, update diagonal part of v[*]
                    if (k > kold) {
                        double temp4 = k * kp1;
                        v(1,k + 1) = 1.0 / temp4;
                        int nsm2 = ns - 2;

                        for (int j = 1; j <= nsm2; j++) {
                            int i = k - j;
                            v(1,i + 1) = v(1,i + 1) - alpha(1,j + 2) * v(1,i + 2);
                        }
                    }

                    // Update V[*] and set W[*]
                    int limit1 = kp1 - ns;
                    double temp5 = alpha(1,ns + 1);

                    for (int iq = 1; iq <= limit1; iq++) {
                        v(1,iq + 1) = v(1,iq + 1) - temp5 * v(1,iq + 2);
                        w(1,iq + 1) = v(1,iq + 1);
                    }

                    g(1,nsp1 + 1) = w(1,2);
                } else {
                    for (int iq = 1; iq <= k; iq++) {
                        double temp3 = iq * (iq + 1);
                        v(1,iq + 1) = 1.0 / temp3;
                        w(1,iq + 1) = v(1,iq + 1);
                    }
                }

                // Compute the g[*] in the work vector w[*]
                int nsp2 = ns + 2;
                if (kp1 >= nsp2) {
                    for (int i = nsp2; i <= kp1; i++) {
                        double limit2 = kp2 - i;
                        double temp0 = alpha(1,i);

                        for (int iq = 1; iq <= limit2; iq++) {
                            w(1,iq + 1) = w(1,iq + 1) - temp0 * w(1,iq + 2);
                        }

                        g(1,i + 1) = w(1,2);
                    }
                }
            }// if K>=NS

            //
            // End block 1
            //

            //
            // Begin block 2
            //
            // Predict a solution p[*], evaluate derivatives using predicted
            // solution, estimate local error at order k and errors at orders
            // k, k-1, k-2 as if constant step size were used.
            //

            // Change phi to phi star

            if (k >= nsp1) {
                for (int i = nsp1; i <= k; i++) {
                    temp1 = beta(1,i + 1);
                    for (int l = 1; l <= n_eqn; l++) {
                        phi(i + 1,l) = temp1 * phi(i + 1,l);
                    }
                }
            }

            // Predict solution and differences
            for (int l = 1; l <= n_eqn; l++) {
                phi(kp2 + 1,l) = phi( kp1 + 1,l);
                phi( kp1 + 1,l) = 0.0;
                p(1, l) = 0.0;
            }

            for (int j = 1; j <= k; j++) {
                int i = kp1 - j;
                int ip1 = i + 1;
                double temp2 = g(1,i + 1);

                for (int l = 1; l <= n_eqn; l++) {
                    p(1,l) = p(1,l) + temp2 * phi(i+1,l);
                    phi(i+1,l) = phi(i+1,l) + phi(ip1+1,l);
                }
            }

            if (nornd) {
                p = y + p * h;
            } else {
                for (int l = 1; l <= n_eqn; l++) {
                    double tau = h * p(1,l) - phi(16,l);
                    p(1,l) = (y)( 1,l) + tau;
                    phi(17,l) = (p(1,l) - (y)(1,l)) - tau;
                }
            }

            double xold = x;
            x = x + h;
            absh = fabs(h);
            yp = func(x, &p);

            // Estimate errors at orders k, k-1, k-2
            erkm2 = 0.0;
            erkm1 = 0.0;
            erk = 0.0;

            for (int l = 1; l <= n_eqn; l++) {
                double temp3 = 1.0 / wt(1,l);
                double temp4 = yp( 1,l) - phi(1 + 1,l);

                if (km2 > 0) {
                    erkm2 = erkm2 + ((phi( km1 + 1,l) + temp4) * temp3) * ((phi( km1 + 1,l) + temp4) * temp3);
                }

                if (km2 >= 0) {
                    erkm1 = erkm1 + ((phi(k + 1,l) + temp4) * temp3) * ((phi( k + 1,l) + temp4) * temp3);
                }

                erk = erk + (temp4 * temp3) * (temp4 * temp3);
            }

            if (km2 > 0) {
                erkm2 = absh * sig(1,km1 + 1) * gstr(1,km2 + 1) * sqrt(erkm2);
            }

            if (km2 >= 0) {
                erkm1 = absh * sig(1,k + 1) * gstr(1,km1 + 1) * sqrt(erkm1);
            }

            double temp5 = absh * sqrt(erk);
            double err = temp5 * (g(1,k + 1) - g(1,kp1 + 1));
            erk = temp5 * sig(1,kp1 + 1) * gstr(1,k + 1);
            knew = k;

            // Test if order should be lowered
            if (km2 > 0) {
                if (fmax(erkm1, erkm2) <= erk) {
                    knew = km1;
                }
            }

            if (km2 == 0) {
                if (erkm1 <= 0.5 * erk) {
                    knew = km1;
                }
            }

            //
            // End block 2
            //

            //
            // If step is successful continue with block 4, otherwise repeat
            // blocks 1 and 2 after executing block 3
            //

            bool success = (err <= epsilon);

            if (!success) {
                //
                // Begin block 3
                //

                // The step is unsuccessful. Restore x, phi[*,*], psi[*]. If
                // 3rd consecutive failure, set order to 1. If step fails more
                // than 3 times, consider an optimal step size. Double error
                // tolerance and return if estimated step size is too small
                // for machine precision.
                //

                // Restore x, phi[*, *] and psi[*]
                phase1 = false;
                x = xold;
                for (int i = 1; i <= k; i++) {
                    temp1 = 1.0 / beta(1,i + 1);
                    int ip1 = i + 1;

                    for (int l = 1; l <= n_eqn; l++) {
                        phi( i + 1,l) = temp1 * (phi(i + 1,l) - phi( ip1 + 1,l));
                    }
                }

                if (k >= 2) {
                    for (int i = 2; i <= k; i++) {
                        psi_(1,i) = psi_(1,i + 1) - h;
                    }
                }

                // On third failure, set order to one.
                // Thereafter, use optimal step size

                ifail = ifail + 1;
                double temp2 = 0.5;
                if (ifail > 3) {
                    if (p5eps < 0.25 * erk) {
                        temp2 = sqrt(p5eps / erk);
                    }
                }

                if (ifail >= 3) {
                    knew = 1;
                }

                h = temp2 * h;
                k = knew;

                if (fabs(h) < fouru * fabs(x)) {
                    crash = true;
                    h = sign(fouru * fabs(x), h);
                    epsilon = epsilon * 2.0;
                    cout << "SALIDA crash3" << endl;
                    return y;                 // Exit
                }

                //
                // End block 3, return to start of block 1
                //
            }

            if (success) {
                break;
            }
        }

        //
        // Begin block 4
        //
        // The step is successful. Correct the predicted solution, evaluate
        // the derivatives using the corrected solution and update the
        // differences. Determine best order and step size for next step.
        //

        kold = k;
        hold = h;

        // Correct and evaluate
        double temp1 = h*g(1,kp1+1);
        if (nornd){
            for (int l = 1; l <= n_eqn; l++) {
                (y)(1,l) = p(1,l) + temp1 * (yp(1,l) - phi(2,l));
            }
        }else{
            for (int l = 1; l <= n_eqn; l++) {
                double rhot =  (yp(1,l) - phi(2,l))*temp1 - phi(17,l);
                (y)(1,l) = p(1,l) + rhot;
                phi(16,l) = ((y)(1,l) - p(1,l)) - rhot;
            }
        }
        yp = func(x,&y);

        // Update differences for next step
        for (int l = 1; l <= n_eqn; l++) {
            phi( kp1 + 1,l) = yp(1,l) - phi( 2,l);
            phi(kp2 + 1,l) = phi( kp1 + 1,l) - phi( kp2 + 1,l);
        }

        for (int i = 1; i <= k; i++) {
            for (int l = 1; l <= n_eqn; l++) {
                phi(i + 1,l) = phi(i + 1,l) + phi(kp1 + 1,l);
            }
        }

        // Estimate error at order k+1 unless
        //- in first phase when always raise order,
        // - already decided to lower order,
        // - step size not constant so estimate unreliable


        double erkp1 = 0.0;

        if ( (knew==km1) || (k!=12) ){
            phase1 = false;
        }

        if (phase1){
            k = kp1;
            erk = erkp1;
        }else{
            if (knew==km1) {
                // lower order
                k = km1;
                erk = erkm1;
            }else{
                if (kp1<=ns){
                    for (int l = 1; l <= n_eqn; l++) {
                        erkp1 = erkp1 + (phi( kp2 + 1,l) / wt(1,l)) * (phi(kp2 + 1,l) / wt(1,l));
                    }

                    erkp1 = gstr(1,kp1+1)*sqrt(erkp1)*absh;
                    // Using estimated error at order k+1, determine
                    // appropriate order for next step

                    if (k>1){
                        if ( erkm1<=fmin(erk,erkp1)) {
                            // lower order
                            k = km1;
                            erk = erkm1;
                        }else{
                            if ( (erkp1<erk) && (k==12) ){
                                // raise order
                                k = kp1;
                                erk = erkp1;
                            }
                        }
                    }else if (erkp1<0.5*erk){
                        // raise order
                        // Here erkp1 < erk < max(erkm1,ermk2) else
                        // order would have been lowered in block 2.
                        // Thus order is to be raised
                        k = kp1;
                        erk = erkp1;
                    }
                } // end if kp1<=ns
            } // end if knew!=km1
        } // end if !phase1

        //With new order determine appropriate step size for next step
        if ( phase1 || (p5eps>=erk*two(1,k+2)) ) {
            hnew = 2.0 * h;
        }else{
            if (p5eps<erk) {
                double temp2 = k + 1;
                double r = p5eps / pow(erk,(1.0 / temp2)); //double r =  pow(p5eps/erk,(1.0 / temp2));
                hnew = absh * fmax(0.5, fmin(0.9, r));
                hnew = sign(fmax(hnew, fouru * fabs(x)), h);
            }else {
                hnew = h;
            }
        }

        h = hnew;

        //
        // End block 4
        //

        // Test for too small tolerances
        if (crash) {
            State_ = DE_STATE::DE_BADACC;
            relerr = epsilon * releps;       // Modify relative and absolute
            abserr = epsilon * abseps;       // accuracy requirements
            y = yy;                   // Copy last step
            t = x;
            told = t;
            OldPermit = true;
            cout << "SALIDA DE_BADACC" << endl;
            return y;                       // Weak failure exit
        }

        nostep = nostep+1;  // Count total number of steps

        // Count number of consecutive steps taken with the order of
        // the method being less or equal to four and test for stiffness

        kle4 = kle4+1;

        if (kold>  4) {
            kle4 = 0;
        }

        if (kle4>=50) {
            stiff = true;
        }
    }
}