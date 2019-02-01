/***************************************************************
* Dynare Program for ECON6201, HW5, Problem 2
* Simulate the Ramsey Model
*
* Paul Gaggl, pgaggl@uncc.edu, 10/30/2015
***************************************************************/

    var y, c, k, i; /*declare endogenous variables*/
    parameters n, g, beta, alpha, sigma, delta; /*declare parameters*/

    /*assigne values to parameters*/
    alpha = 1/3;
    delta = 1;
    g=0.02;
    n=0.01;
    sigma=1;
    beta=1/1.01;

/***************************************************************
* Equilibrium Conditions for the Model
* Hint: Remember that, for dynare, capital needs to be
* defined as "end of period". That is,
*      k  = captial at the end of t
*   k(-1) = captial at the end of t-1 (=beginning of t)
***************************************************************/
    
    model;
        y = k(-1)^alpha;
        ((1+n)*(1+g)^sigma)*c^(-sigma) = (1-delta+alpha*(y(+1)/k))*beta*c(+1)^(-sigma);
        i = y-c;
        (1+n)*(1+g)*k = (1-delta)*k(-1) + i;
    end;

/***************************************************************
* Initial values. Here I use the theoretical steady state
* as initial values.
***************************************************************/
    
    initval;
        k=0.5*((alpha*beta)/((1+n)*(1+g)^sigma - beta*(1-delta)))^(1/(1-alpha));
        y=k^alpha;
        i = ((1+n)*(1+g) - (1-delta))*k;
        c=y-i;
    end;
    
    check; //check whether the system is stable.                         
           //If this gives you an error message this means that          
           //something is wrong with your equations in the "model"       
           //block.    

    
/***************************************************************
* New Situation after unanticipated permanent shock
***************************************************************/

    endval;
        k=((alpha*beta)/((1+n)*(1+g)^sigma - beta*(1-delta)))^(1/(1-alpha));
        y=k^alpha;
        i = ((1+n)*(1+g) - (1-delta))*k;
        c=y-i;
    end;
    steady; //compute new steady state
    

/***************************************************************
* Deterministic Simulation
***************************************************************/
    
    //this command simulated the model for 100 periods          
    //starting with the initial values in the "initval" block,  
    //and simulating the model with the new, differnt z,        
    //given in the "endval" block.                              
    simul(
        periods=500
    );
    
    
