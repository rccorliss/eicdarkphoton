//                                                                    -*-C++-*-

/*
Calculating the A' production cross section on heavy target

(1)             (2)     	   
      	 .                     .
      	. A'                  . A'
       .                     .	  
  e   .	      e'    e	    .	e'
  ---o---o----->    ---o---o----->	   
        (   	      (  	  	   
         )	       )          
  ------o------>    --o---------->
  Z           Z'    Z           Z'
*/

// Parameter:
//   e_in:      incoming electron e
//   e_out:     outgoing electron e'
//   q_out:     outgoing A'
//   mA:        mass of A'

long double DMHeavyCS(const FourVector &e_in,  const FourVector &e_out,
		 const FourVector &q_out, const long double mA);

/*
Calculating the four QED Background graphs:

(1)   |e+         (2)      |e+       (3)	       (4)		   
      |	  e-	           |  e-       -----o----->      ------o----->
      o--->	           o----       e     )   e'	 e      )   e'
     (	                   )	            (    e+	       o   /e+
  e   )	      e'    e	  (	e'           o-----            |\ /
  ---o---o----->    ---o---o----->           | 		       | X	   
        (   	      (  	             o-----	       |/ \	   
         )	       )                    (    e-	       o   \e-
  ------o------>    --o---------->     Z     )	 Z'      Z    (	    Z'
  Z           Z'    Z           Z'     -----o----->      ------o----->

*/

// Parameter:
//   e_in:      incoming electron e
//   e_out:     outgoing electron e'
//   q_out:     outgoing A', i.e. the sum of the leptons: A' = e+ + e'
//   mA:        mass of A', i.e. sqrt((e+ + e')^2)
//   theta_e12:
//   phi_e12:   A' decay angles in A' rest frame

long double QEDBackground(const FourVector &e_in,  const FourVector &e_out,
		     const FourVector &q_out, const long double mA, 
		     const long double theta_e12, const long double phi_e12);
