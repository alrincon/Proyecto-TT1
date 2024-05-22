#include "../include/DEInteg.h"

Matrix extractRowt(Matrix* input, int column){
    Matrix res(input->getFilas(), 1);
    for(int i = 1; i <= input->getFilas(); i++){
        res(i, 1) = (*input)(i, column);
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

Matrix DEInteg(Matrix (*func)(double, Matrix*), double t, double tout, double relerr, double abserr, int n_eqn, Matrix *y) {
    cout << "entra" << endl;

    double twou = 2 * 2.2204e-16;
    double fouru = 4 * 2.2204e-16;

    int State_ = DE_STATE::DE_INIT;
    bool PermitTOUT = true; // Allow integration past tout by default
    bool OldPermit;
    double told = 0;

    Matrix two(14, 1);
    for (int i = 1; i <= 14; i++) {
        two(i, 1) = pow(2.0, i-1);
    }

    cout << "checkpoint 1" << endl;
    Matrix gstr(14, 1);
    gstr(1, 1) = 1.0;
    gstr(2, 1) = 0.5;
    gstr(3, 1) = 0.0833;
    gstr(4, 1) = 0.0417;
    gstr(5, 1) = 0.0264;
    gstr(6, 1) = 0.0188;
    gstr(7, 1) = 0.0143;
    gstr(8, 1) = 0.0114;
    gstr(9, 1) = 0.00936;
    gstr(10, 1) = 0.00789;
    gstr(11, 1) = 0.00679;
    gstr(12, 1) = 0.00592;
    gstr(13, 1) = 0.00524;
    gstr(14, 1) = 0.00468;

    //Hay que darle la vuelta a todos
    Matrix yy(n_eqn, 1);    // Allocate vectors with proper dimension
    Matrix wt(n_eqn, 1);
    Matrix p(n_eqn, 1);
    Matrix yp(n_eqn, 1);
    Matrix phi(n_eqn, 17);
    Matrix g(14, 1);
    Matrix sig(14, 1);
    Matrix rho(14, 1);
    Matrix w(13, 1);
    Matrix alpha(13, 1);
    Matrix beta(13, 1);
    Matrix v(13, 1);
    Matrix psi_(13, 1);

    if (t == tout) return *y;

    cout << "checkpoint 2" << endl;
    double epsilon = std::max(relerr, abserr);
    if (relerr < 0.0 || abserr < 0.0 || epsilon <= 0.0 || State_ > DE_STATE::DE_INVPARAM || (State_ != DE_STATE::DE_INIT && t != told)) {
        State_ = DE_STATE::DE_INVPARAM;
        return *y;
    }

    double del = tout - t;
    double absdel = abs(del);

    double tend = t + 100.0 * del;
    if (!PermitTOUT) tend = tout;

    int nostep = 0;
    int kle4 = 0;
    bool stiff = false;
    double releps = relerr / epsilon;
    double abseps = abserr / epsilon;

    bool start = true;
    double x = t;
    double delsgn = sign(1.0, del);
    double h = sign(max(fouru * abs(x), abs(tout - x)), tout - x);

    cout << "checkpoint 3" << endl;
    if  ( (State_==DE_STATE::DE_INIT) || (delsgn*del<=0.0) ){
        // On start and restart also set the work variables x and yy(*),
        // store the direction of integration and initialize the step size
        cout << "entra" << endl;
        start  = true;
        x      = t;
        cout << "entra" << endl;
        yy     = (*y);
        cout << "entra" << endl;
        delsgn = sign(1.0, del);
        cout << "entra" << endl;
        h      = sign( max(fouru*abs(x), abs(tout-x)), tout-x );
    }


    cout << "checkpoint 4" << endl;
    while (true) {
        Matrix yout(n_eqn, 1);
        Matrix ypout(n_eqn, 1);
        double hi;
        int ki;
        int kold = 0;

        // If already past output point, interpolate solution and return
        if (abs(x - t) >= absdel) {
            g(2, 1) = 1.0;
            rho(2, 1) = 1.0;
            hi = tout - x;
            ki = kold + 1;

            // Initialize w[*] for computing g[*]
            for (int i = 1; i <= ki; i++) {
                double temp1 = i;
                w(i + 1, 1) = 1.0 / temp1;
            }

            // Compute g[*]
            double term = 0.0;

            for (int j = 2; j <= ki; j++) {
                double psijm1 = psi_(j, 1);
                double gamma = (hi + term) / psijm1;
                double eta = hi / psijm1;

                for (int i = 1; i <= ki + 1 - j; i++) {
                    w(i + 1, 1) = gamma * w(i + 1, 1) - eta * w(i + 2, 1);
                }

                g(j + 1, 1) = w(2, 1);
                rho(j + 1, 1) = gamma * rho(j, 1);
                term = psijm1;
            }


            cout << "checkpoint 5" << endl;
            // Interpolate for the solution yout and for
            // the derivative of the solution ypout
            for (int j = 1; j <= ki; j++) {
                int i = ki + 1 - j;
                yout = yout + extractRowt(&phi, i + 1) * g(i + 1, 1);
                ypout = ypout + extractRowt(&phi, i + 1) * rho(i + 1, 1);
            }

            yout = *y + yout * hi;
            y = &yout;

            State_ = DE_STATE::DE_DONE; // Set return code
            t = tout;             // Set independent variable
            told = t;                // Store independent variable
            OldPermit = PermitTOUT;
            return *y;                       // Normal exit
        }

        cout << "checkpoint 6" << endl;
        // If cannot go past output point and sufficiently close,
        // extrapolate and return
        if (!PermitTOUT && (abs(tout - x) < fouru * abs(x))) {
            h = tout - x;
            yp = func(x, &yy);          // Compute derivative yp(x)
            Matrix typ(n_eqn,1);
            typ = yy + (yp) * h;
            y = &typ;                // Extrapolate vector from x to tout
            State_ = DE_STATE::DE_DONE; // Set return code
            t = tout;             // Set independent variable
            told = t;                // Store independent variable
            OldPermit = PermitTOUT;
            return *y;                       // Normal exit
        }

        // Limit step size, set weight vector and take a step
        h = sign(min(abs(h), abs(tend - x)), h);
        for (int l = 1; l <= n_eqn; l++) {
            wt(l, 1) = releps * abs(yy(l, 1)) + abseps;
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

        double crash;

        cout << "checkpoint 7" << endl;
        if (abs(h) < fouru * abs(x)) {
            h = sign(fouru * abs(x), h);
            crash = true;
            return *y;           // Exit
        }

        double p5eps = 0.5 * epsilon;
        crash = false;
        g(2, 1) = 1.0;
        g(3, 1) = 0.5;
        sig(2, 1) = 1.0;

        double ifail = 0;

        // If error tolerance is too small, increase it to an
        // acceptable value.

        double round = 0.0;

        for (int l = 1; l <= n_eqn; l++) {
            round = round + ((*y)(l, 1) * (*y)(l, 1)) / (wt(l, 1) * wt(l, 1));
        }

        round = twou * sqrt(round);

        if (p5eps < round) {
            epsilon = 2.0 * round * (1.0 + fouru);
            crash = true;
            return *y;
        }


        double absh;
        int k;
        double hold;
        double hnew;
        bool phase1;
        bool nornd;

        cout << "checkpoint 8" << endl;
        if (start) {
            // Initialize. Compute appropriate step size for first step.
            yp = func(x, y);
            double sum = 0.0;

            cout << "checkpoint 81" << endl;
            for (int l = 1; l <= n_eqn; l++) {
                phi(l, 2) = yp(l, 1);
                phi(l, 3) = 0.0;
                sum = sum + (yp(l, 1) * yp(l, 1)) / (wt(l, 1) * wt(l, 1));

                cout << "checkpoint 81" << endl;
            }

            sum = sqrt(sum);
            absh = abs(h);

            if (epsilon < 16.0 * sum * h * h) {
                absh = 0.25 * sqrt(epsilon / sum);
            }

            cout << "checkpoint 82" << endl;
            h = sign(max(absh, fouru * abs(x)), h);
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
                    phi(l, 16) = 0.0;
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

        cout << "checkpoint 9" << endl;
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
                beta(ns + 1, 1) = 1.0;
                realns = ns;
                alpha(ns + 1, 1) = 1.0 / realns;
                temp1 = h * realns;
                sig(nsp1 + 1, 1) = 1.0;

                if (k >= nsp1) {
                    for (int i = nsp1; i <= k; i++) {
                        double im1 = i - 1;
                        double temp2 = psi_(im1 + 1, 1);
                        psi_(im1 + 1, 1) = temp1;
                        beta(i + 1, 1) = beta(im1 + 1, 1) * psi_(im1 + 1, 1) / temp2;
                        temp1 = temp2 + h;
                        alpha(i + 1, 1) = h / temp1;
                        int reali = i;
                        sig(i + 2, 1) = reali * alpha(i + 1, 1) * sig(i + 1, 1);
                    }
                }

                psi_(k + 1, 1) = temp1;

                // Compute coefficients g[*]; initialize v[*] and set w[*].
                if (ns > 1) {
                    // If order was raised, update diagonal part of v[*]
                    if (k > kold) {
                        double temp4 = k * kp1;
                        v(k + 1, 1) = 1.0 / temp4;
                        int nsm2 = ns - 2;

                        for (int j = 1; j <= nsm2; j++) {
                            int i = k - j;
                            v(i + 1, 1) = v(i + 1, 1) - alpha(j + 2, 1) * v(i + 2, 1);
                        }
                    }

                    // Update V[*] and set W[*]
                    int limit1 = kp1 - ns;
                    double temp5 = alpha(ns + 1, 1);

                    for (int iq = 1; iq <= limit1; iq++) {
                        v(iq + 1, 1) = v(iq + 1, 1) - temp5 * v(iq + 2, 1);
                        w(iq + 1, 1) = v(iq + 1, 1);
                    }

                    g(nsp1 + 1, 1) = w(2, 1);
                } else {
                    for (int iq = 1; iq <= k; iq++) {
                        double temp3 = iq * (iq + 1);
                        v(iq + 1, 1) = 1.0 / temp3;
                        w(iq + 1, 1) = v(iq + 1, 1);
                    }
                }

                // Compute the g[*] in the work vector w[*]
                int nsp2 = ns + 2;
                if (kp1 >= nsp2) {
                    for (int i = nsp2; i <= kp1; i++) {
                        double limit2 = kp2 - i;
                        double temp0 = alpha(i, 1);

                        for (int iq = 1; iq <= limit2; iq++) {
                            w(iq + 1, 1) = w(iq + 1, 1) - temp0 * w(iq + 2, 1);
                        }

                        g(i + 1, 1) = w(2, 1);
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
                    temp1 = beta(i + 1, 1);
                    for (int l = 1; l <= n_eqn; l++) {
                        phi(l, i + 1) = temp1 * phi(l, i + 1);
                    }
                }
            }

            // Predict solution and differences
            for (int l = 1; l <= n_eqn; l++) {
                phi(l, kp2 + 1) = phi(l, kp1 + 1);
                phi(l, kp1 + 1) = 0.0;
                p(l, 1) = 0.0;
            }

            for (int j = 1; j <= k; j++) {
                int i = kp1 - j;
                int ip1 = i + 1;
                double temp2 = g(i + 1, 1);

                for (int l = 1; l <= n_eqn; l++) {
                    p(l, 1) = p(l, 1) + temp2 * phi(l, i + 1);
                    phi(l, i + 1) = phi(l, i + 1) + phi(l, ip1 + 1);
                }
            }

            if (nornd) {
                p = *y + p * h;
            } else {
                for (int l = 1; l <= n_eqn; l++) {
                    double tau = h * p(l, 1) - phi(l, 16);
                    p(l, 1) = (*y)(l, 1) + tau;
                    phi(l, 17) = (p(l, 1) - (*y)(l, 1)) - tau;
                }
            }

            double xold = x;
            x = x + h;
            absh = abs(h);
            yp = func(x, &p);

            // Estimate errors at orders k, k-1, k-2
            erkm2 = 0.0;
            erkm1 = 0.0;
            erk = 0.0;

            for (int l = 1; l <= n_eqn; l++) {
                double temp3 = 1.0 / wt(l, 1);
                double temp4 = yp(l, 1) - phi(l, 1 + 1);

                if (km2 > 0) {
                    erkm2 = erkm2 + ((phi(l, km1 + 1) + temp4) * temp3) * ((phi(l, km1 + 1) + temp4) * temp3);
                }

                if (km2 >= 0) {
                    erkm1 = erkm1 + ((phi(l, k + 1) + temp4) * temp3) * ((phi(l, k + 1) + temp4) * temp3);
                }

                erk = erk + (temp4 * temp3) * (temp4 * temp3);
            }

            if (km2 > 0) {
                erkm2 = absh * sig(km1 + 1, 1) * gstr(km2 + 1, 1) * sqrt(erkm2);
            }

            if (km2 >= 0) {
                erkm1 = absh * sig(k + 1, 1) * gstr(km1 + 1, 1) * sqrt(erkm1);
            }

            double temp5 = absh * sqrt(erk);
            double err = temp5 * (g(k + 1, 1) - g(kp1 + 1, 1));
            erk = temp5 * sig(kp1 + 1, 1) * gstr(k + 1, 1);
            knew = k;

            // Test if order should be lowered
            if (km2 > 0) {
                if (max(erkm1, erkm2) <= erk) {
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
                    temp1 = 1.0 / beta(i + 1, 1);
                    int ip1 = i + 1;

                    for (int l = 1; l <= n_eqn; l++) {
                        phi(l, i + 1) = temp1 * (phi(l, i + 1) - phi(l, ip1 + 1));
                    }
                }

                if (k >= 2) {
                    for (int i = 2; i <= k; i++) {
                        psi_(i, 1) = psi_(i + 1, 1) - h;
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

                if (abs(h) < fouru * abs(x)) {
                    crash = true;
                    h = sign(fouru * abs(x), h);
                    epsilon = epsilon * 2.0;
                    return *y;                 // Exit
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
        double temp1 = h*g(kp1+1,1);
        if (nornd){
            for (int l = 1; l <= n_eqn; l++) {
                (*y)(l,1) = p(l,1) + temp1 * (yp(l,1) - phi(l, 2));
            }
        }else{
            for (int l = 1; l <= n_eqn; l++) {
                double rhot =  (yp(l,1) - phi(l, 2))*temp1 - phi(l, 17);
                (*y)(l,1) = p(l,1) + rhot;
                phi(l, 16) = ((*y)(l,1) - p(l,1)) - rhot;
            }
        }

        yp = func(x,y);

        // Update differences for next step
        for (int l = 1; l <= n_eqn; l++) {
            phi(l, kp1 + 1) = yp(l,1) - phi(l, 2);
            phi(l, kp2 + 1) = phi(l, kp1 + 1) - phi(l, kp2 + 1);
        }

        for (int i = 1; i <= k; i++) {
            for (int l = 1; l <= n_eqn; l++) {
                phi(l, i + 1) = phi(l, i + 1) + phi(l, kp1 + 1);
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
                        erkp1 = erkp1 + (phi(l, kp2 + 1) / wt(l,1)) * (phi(l, kp2 + 1) / wt(l,1));
                    }

                    erkp1 = gstr(kp1+1,1)*sqrt(erkp1)*absh;
                    // Using estimated error at order k+1, determine
                    // appropriate order for next step

                    if (k>1){
                        if ( erkm1<=min(erk,erkp1)) {
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
        if ( phase1 || (p5eps>=erk*two(k+2,1)) ) {
            hnew = 2.0 * h;
        }else{
            if (p5eps<erk) {
                double temp2 = k + 1;
                double r = p5eps / pow(erk,(1.0 / temp2)); //double r =  pow(p5eps/erk,(1.0 / temp2));
                hnew = absh * max(0.5, min(0.9, r));
                hnew = sign(max(hnew, fouru * abs(x)), h);
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
            y = &yy;                   // Copy last step
            t = x;
            told = t;
            OldPermit = true;
            return *y;                       // Weak failure exit
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