$title A-AUGMECON2:Advanced version of AUGMECON2 for Multiobjective Optimization
$ontext
The advanced version of AUGMECON2, denoted as A-AUGMECON2, is provided with
a novel pruning algorithm that avoids solving redundant optimizations,
leading to a significant reduction of the computational burden.

The algorithm has been developed starting from the AUGMECON2 code
available in GAMS Model Library at the following link:
https://www.gams.com/latest/gamslib_ml/libhtml/gamslib_epscmmip.html

For the sake of simplicity, A-AUGMECON2 is here applied to a trivial
Multi-Objective Integer Programming problem(specifically a Multi-Objective
Multi-Dimensional Knapsack Problem) with 50 binary variables X, 2 objective
functions and 2 constraints. The higher the complexity of the problem to be
optimized in terms of number of objective functions and density of the Pareto
frontier, the more the novel A-AUGMECON2 outperforms the standard AUGMECON2.


Additional information can be found at:
https://doi.org/10.1016/j.apenergy.2021.117283

Marina Petrelli, Davide Fioriti, Alberto Berizzi, Cristian Bovo, Davide Poli,
"A novel multi-objective method with online Pareto pruning for multi-year
optimization of rural microgrids", Applied Energy, Volume 299, 2021, 117283.


The paper compares the results of the novel A-AUGMECON2 with AUGMECON2 applied
to a rural microgrid planning problem, evaluating costs, emissions, land
use, job creation and public lighting coverage. A-AUGMECON2 allows to reduce
the computational burden by 48% and the number of points in the Pareto frontier
by 42%, removing redundant simulations and improving the quality and readibility
of the results.


INSTRUCTIONS FOR REPLICATION:
In order to apply A-AUGMECON2 to a different model, it is enough to substitute
the Knapsack Problem (whose description finishes at line 98) with the desired
optimization problem. It is necessary to specify the set of objective functions
K, the related variable z(K) and direction of optimization dir(K). The desired
density of the Pareto curve should be set at line 198.
$offtext


$eolcom //
$STitle Example model definitions

Sets
   I       constraints         /  i1* i2 /
   J       decision variables  /  j1*j50 /
   K       objective functions /  k1* k2 /

$set min -1
$set max +1
Parameter
   dir(K) direction of the objective functions / k1 %max%, k2 %max% /
   b(I)    RHS of the constraints / i1 1445, i2 1502.5 /

Table c(K,J) matrix of objective function coefficients C
    j1  j2  j3  j4  j5  j6  j7  j8  j9 j10 j11 j12 j13 j14 j15 j16 j17
k1  21  69  26  92  77  30  96  80  60  61  52  92  19  10  63  34 100
k2  24  92  53  25  10  31  83  34  64  69  95  40  59  87  13  94  53
+
   j18 j19 j20 j21 j22 j23 j24 j25 j26 j27 j28 j29 j30 j31 j32 j33 j34
k1  60  11  12  37 100  74  17  60  69  49  69  49  59  17  21  74  85
k2  52  61  53  78  34  89  32  28  56  52  40  41  59  35  96  72  55
+
   j35 j36 j37 j38 j39 j40 j41 j42 j43 j44 j45 j46 j47 j48 j49 j50
k1  83  41  29  63  56  38  66  92  25  84  89  21  46  94  96  92
k2 100  44  90  66  59  22  72  25  36  16  56  91  61  56  66  53
;

Table a(I,J) matrix of technological coefficients A
    j1  j2  j3  j4  j5  j6  j7  j8  j9 j10 j11 j12 j13 j14 j15 j16 j17
i1  84  49  68  20  97  74  60  30  13  95  19  41  17  95  73  12  66
i2  19  96  93  64  72  91  32  96  44  76  69  82  51  38  52  22  83
+
   j18 j19 j20 j21 j22 j23 j24 j25 j26 j27 j28 j29 j30 j31 j32 j33 j34
i1  55  75  20  56  80  59  66  25  70  95  96  62  74  31  59  21  85
i2  27  70  56  29  89  86  48  13  95  66  94  16  44  67  90  48  29
+
   j35 j36 j37 j38 j39 j40 j41 j42 j43 j44 j45 j46 j47 j48 j49 j50
i1  45  97  23  53  51  95  58  68  62  45  83  82  47  15  52  72
i2  90  54  77  28 100  86  51  62  40  54  21  55  50  62  51  77
;

Variables
   Z(K)      objective function variables
   X(J)      decision variables
Binary Variables X;

Equations
   objfun(K) objective functions
   con(I)    constraints;

objfun(K)..   sum(J, c(K,J)*X(J)) =e= Z(K);
con(I)..      sum(J, a(I,J)*X(J)) =l= b(I);

Model example / all /;

$STitle eps-constraint method

Set k1(k)  the first element of k
    km1(k) all but the first elements of k
    kk(k)  active objective function in constraint allobj;

k1(k)$(ord(k)=1) = yes; km1(k)=yes; km1(k1) = no;

Parameter
    rhs(k)    right hand side of the constrained obj functions in eps-constraint
    maxobj(k) maximum value from the payoff table
    minobj(k) minimum value from the payoff table
    numk(k)   ordinal value of k starting with 1
    val(k)    desired value of objective function in payoff table
    price(k)  penalty for difference from desired value in payoff table;

val(k)=0; price(k)=0;

Scalar
   iter         total number of iterations
   elapsed_time elapsed time for eps-constraint
   start        start time
   finish       finish time

Variables
   a_objval   auxiliary variable for the objective function
   obj        auxiliary variable during the construction of the payoff table

Positive Variables
   s(k)       slack or surplus variables for the eps-constraints
   s1(k)      slack or surplus variables for the payoff table
   s2(k)      slack or surplus variables for the payoff table

Equations
   allobj     all the objective functions in one expression
   con_payoff(k) constraint on desired value of objective function in payoff table
   con_obj(k) constrained objective functions
   augm_obj   augmented objective function to avoid weakly efficient solutions;

con_obj(km1)..   z(km1) - dir(km1)*s(km1) =e= rhs(km1);

* We optimize the first objective function and put the others as constraints
* the second term is for avoiding weakly efficient points

augm_obj.. a_objval =e= sum(k1,dir(k1)*z(k1)/(maxobj(k1)-minobj(k1)))
    + 1e-3*sum(km1,power(10,-(numk(km1)-1))*s(km1)/(maxobj(km1)-minobj(km1)));

allobj.. sum(kk, dir(kk) * z(kk)) - sum(k, price(k) * (s1(k) + s2(k))) =e= obj ;

con_payoff(k).. z(k) + s1(k) - s2(k) =E= val(k) ;

Model mod_payoff    / example, allobj, con_payoff / ;
Model mod_epsmethod / example, con_obj, augm_obj / ;

Parameter
   payoff(k,k)  payoff tables entries;
Alias(k,kp);
Alias(k,k2);

option optcr=0, limrow=0, limcol=0, solprint=off, solvelink=%Solvelink.LoadLibrary%;

* Generate payoff table applying lexicographic optimization
loop(kp,
* Modify the priority order to exploit these results as Pareto points, removing redundant simulations
  loop(k2,
    if (ord(k2) = 1,
      kk(kp) = yes;
    elseif ord(k2) <= ord(kp),
      kk(k2-1) = yes;
    else
      kk(k2) = yes;
    );
    solve mod_payoff using mip maximizing obj;
    val(kk) = z.l(kk); // desired value of the objective function
    price(kk) = 1e5 ; // big constant to penalize distance from desired value
    kk(k) = no;
  );
  payoff(kp,k) = z.l(k);
  kk(k) = no;
* release the desired values of the objective functions for the new iteration
  val(k) = 0; price(k) = 0;
);


if (mod_payoff.modelstat<>%ModelStat.Optimal% and
    mod_payoff.modelstat<>%ModelStat.Integer Solution%,
   abort 'no optimal solution for mod_payoff');

file fx  / 2kp50_augmecon2_results.txt /;
put fx ' PAYOFF TABLE'/   ;
loop (kp,
   loop(k, put payoff(kp,k):12:2);
   put /);

minobj(k)=smin(kp,payoff(kp,k));
maxobj(k)=smax(kp,payoff(kp,k));


$if not set gridpoints $set gridpoints 491
*$if not set gridpointsred $set gridpointsred 200
* gridpointsred coould be used to set reduced densities for some objectives
Set g            grid points /g0*g%gridpoints%/
*    gr_red(gr) reduced grid points /gr0*gr%gridpointsred%/
    grid(k,g)    grid
Parameter
    gridrhs(k,g) rhs of eps-constraint at grid point
    maxg(k)      maximum point in grid for objective
    posg(k)      grid position of objective
    firstOffMax, lastZero, payoffSol some counters
    numg(g)      ordinal value of g starting with 0
    step(k)      step of grid points in objective functions
    jump(k)      jumps in the grid points traversing;

lastZero=1; loop(km1, numk(km1)=lastZero; lastZero=lastZero+1); numg(g) = ord(g)-1;

* Here we could define different grid intervals for different objectives
* using set gr_red
grid(km1,g) = yes;
maxg(km1)   = smax(grid(km1,g), numg(g));
step(km1)   = (maxobj(km1)- minobj(km1))/maxg(km1);
gridrhs(grid(km1,g))$(dir(km1)=-1) = maxobj(km1) - numg(g)/maxg(km1)*(maxobj(km1)- minobj(km1));
gridrhs(grid(km1,g))$(dir(km1)= 1) = minobj(km1) + numg(g)/maxg(km1)*(maxobj(km1)- minobj(km1));

put / ' Grid points' /;
loop (g,
   loop(km1, put gridrhs(km1,g):12:2);
   put /);
put / 'Efficient solutions' /;

* Walk the grid points and take shortcuts if the model becomes infeasible or
* if the calculated slack variables are greater than the step size

$eval vNTOT power(%gridpoints% + 1, card(K)-1)

set        num     /1*%vNTOT%/
           n(num) ;
parameter  numn(num)  ordinal value of num starting with 1
           numord     ordinal number of the gridpoint examined
           v(num)     vector of points to be examined; //v(n)=1 if point to be solved, v(n)=0 if point to be skipped
parameter  Nc         number of combinations to be excluded from examination because redundant
           comb       counter of combinations
           combb      scalar for computing combinations
           deltav(k)  identifier of solutions to be excluded
           posg_v(k)  position to be excluded
           numordp(k) ordinal number of payoff table points
           posgp(k)   grid positions of payoff table points
           numordPre
           numordPost;

alias(km1,kmm1); alias (km1,kkm1);

n(num) = yes;
lastZero=1; loop(n, numn(n)=lastZero; lastZero=lastZero+1);
v(n)=1; // initialize all as 1

* this loop computes gridpoints corresponding to payoff table points
loop(km1,  posg(km1)=maxg(km1); posg(kmm1)$(numk(kmm1) NE numk(km1))=0;
           numordp(km1)= sum(kkm1$(numk(kkm1)>1), posg(kkm1)* prod(kmm1$(numk(kmm1) LE numk(kkm1)-1), maxg(kmm1)+1 ))
           + sum(kkm1$(numk(kkm1)=1), posg(kkm1)) + 1;
     );

posg(km1) = 0; iter=0; start=jnow;

repeat

  numord = sum(km1$(numk(km1)>1), posg(km1)* prod(kmm1$(numk(kmm1) LE numk(km1)-1), maxg(kmm1)+1 ))
           + sum(km1$(numk(km1)=1), posg(km1)) + 1;

  payoffSol=0;

  if(sum(n$(numn(n)=numord), v(n)) = 1,

* the following two loops verify if the current gridpoint corresponds to a point of the payofftable
* if so, the optimization is skipped and the solution of the payofftable is directly included
         loop(k1$(sum(km1, posg(km1)) = 0),
                z.l(k) = payoff(k1,k);
                payoffSol=1;
                iter=iter+1;
                put iter:5:0;
                loop(k, put z.l(k):12:2);
                put /;
                jump(km1)$(%min%=dir(km1))=1+floor( (maxobj(km1) - z.l(km1)) / step(km1) ); // posg(km1)=0
                jump(km1)$(%max%=dir(km1))=1+floor( (z.l(km1) - minobj(km1)) / step(km1) );
              );

         loop(km1$(posg(km1)=maxg(km1) and sum(kmm1$(numk(kmm1) NE numk(km1)), posg(kmm1))=0),
                z.l(k) = payoff(km1,k);
                payoffSol=1;
                iter=iter+1;
                put iter:5:0;
                loop(k, put z.l(k):12:2);
                put /;
                jump(km1)$(%min%=dir(km1))=1+floor( (minobj(km1) - z.l(km1)) / step(km1) ); // posg(km1)=maxg(km1)
                jump(km1)$(%max%=dir(km1))=1+floor( (z.l(km1) - maxobj(km1)) / step(km1) );
                jump(kmm1)$(%min%=dir(kmm1) and (numk(kmm1) NE numk(km1)))=1+floor( (maxobj(kmm1) - z.l(kmm1)) / step(kmm1) ); // posg(kmm1)=0
                jump(kmm1)$(%max%=dir(kmm1) and (numk(kmm1) NE numk(km1)))=1+floor( (z.l(kmm1) - minobj(kmm1)) / step(kmm1) );
              );

     if(payoffSol=1,
         Nc = prod(km1, jump(km1));
         for( comb = 0 to Nc-1 by 1, // find the positions of the points to be skipped
              combb = comb ;
              loop(km1, deltav(km1) = mod(combb,jump(km1));
                        combb = floor(combb/jump(km1));
                        posg_v(km1) = posg(km1) + deltav(km1);
                   );
               numord = sum(km1$(numk(km1)>1), posg_v(km1)* prod(kmm1$(numk(kmm1) LE numk(km1)-1), maxg(kmm1)+1 ))
                        + sum(km1$(numk(km1)=1), posg_v(km1)) + 1;
               v(n)$(numn(n)=numord) = 0 ;
             );
     else
         rhs(km1) = sum(grid(km1,g)$(numg(g)=posg(km1)), gridrhs(km1,g));
         solve mod_epsmethod maximizing a_objval using mip;
         if (mod_epsmethod.modelstat<>%ModelStat.Optimal% and
             mod_epsmethod.modelstat<>%ModelStat.Integer Solution%,
             lastZero = 0; loop(km1$(posg(km1)>0 and lastZero=0), lastZero=numk(km1));
             numordPre = sum(km1$(numk(km1)>1), posg(km1)* prod(kmm1$(numk(kmm1) LE numk(km1)-1), maxg(kmm1)+1 ))
                         + sum(km1$(numk(km1)=1), posg(km1)) + 1;
             posg(km1)$(numk(km1)<=lastZero) = maxg(km1); // skip all solves for more demanding values of rhs(km1)
             numordPost = sum(km1$(numk(km1)>1), posg(km1)* prod(kmm1$(numk(kmm1) LE numk(km1)-1), maxg(kmm1)+1 ))
                         + sum(km1$(numk(km1)=1), posg(km1)) + 1;
             loop(km1$((numordPre<numordp(km1))and (numordPost>=numordp(km1))), // include points of the payofftable even if the counter would go beyond
                  z.l(k) = payoff(km1,k);
                  payoffSol=1;
                  iter=iter+1;
                  put iter:5:0 ,';';
                  loop(k, put z.l(k):18:5 ,';');
                  put /;
                  jump(km1)$(%min%=dir(km1))=1+floor( (minobj(km1) - z.l(km1)) / step(km1) ); // posg(km1)=maxg(km1)
                  jump(km1)$(%max%=dir(km1))=1+floor( (z.l(km1) - maxobj(km1)) / step(km1) );
                  jump(kmm1)$(%min%=dir(kmm1) and (numk(kmm1) NE numk(km1)))=1+floor( ( maxobj(kmm1) - z.l(kmm1)) / step(kmm1) ); // posg(kmm1)=0
                  jump(kmm1)$(%max%=dir(kmm1) and (numk(kmm1) NE numk(km1)))=1+floor( (z.l(kmm1) - minobj(kmm1)) / step(kmm1) );
                 );
                 if(payoffSol=1,
                  Nc = prod(km1, jump(km1)); display jump, Nc;
                  for( comb = 0 to Nc-1 by 1, // find the positions of the points to be skipped
                       combb = comb ;
                       loop(km1, deltav(km1) = mod(combb,jump(km1));
                                 combb = floor(combb/jump(km1));
                                 posg_v(km1) = posg(km1) + deltav(km1);
                            );
                        numord = sum(km1$(numk(km1)>1), posg_v(km1)* prod(kmm1$(numk(kmm1) LE numk(km1)-1), maxg(kmm1)+1 ))
                                 + sum(km1$(numk(km1)=1), posg_v(km1)) + 1;
                        v(n)$(numn(n)=numord) = 0 ;
                      );
                    );
         else
             iter=iter+1;
             put iter:5:0;
             loop(k, put z.l(k):12:2);
             put /;
             jump(km1)=1+floor(s.L(km1)/step(km1));
             Nc = prod(km1, jump(km1));
             for ( comb = 0 to Nc-1 by 1, // find the positions of the points to be skipped
                   combb = comb ;
                   loop(km1, deltav(km1) = mod(combb,jump(km1));
                             combb = floor(combb/jump(km1));
                             posg_v(km1) = posg(km1) + deltav(km1);
                         );
                   numord = sum(km1$(numk(km1)>1), posg_v(km1)* prod(kmm1$(numk(kmm1) LE numk(km1)-1), maxg(kmm1)+1 ))
                            + sum(km1$(numk(km1)=1), posg_v(km1)) + 1;
                   v(n)$(numn(n)=numord) = 0 ;
                 );
             );  // end of if on modelstat

      ); // end of if on payoffSol

   );  // end of if on v(n)


* Proceed forward in the grid
  firstOffMax = 0;
  loop(km1$(posg(km1)<maxg(km1) and firstOffMax=0),
     posg(km1)=posg(km1)+1; firstOffMax=numk(km1));
  posg(km1)$(numk(km1)<firstOffMax) = 0;
until sum(km1$(posg(km1)=maxg(km1)),1)= card(km1) and firstOffMax=0;


finish=jnow; elapsed_time=(finish-start)*60*60*24;

put /;
put 'Elapsed time: ',elapsed_time:10:2, ' seconds' / ;
