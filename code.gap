###
##  This code for GAP ( version 4.14.0 of 2024-12-05 ) accompanies the paper
##  
##  "On generation by triality automorphisms"
## 
##  by Danila O. Revin and Andrei V. Zavarnitsine
##
##  Date: July, 2025

##  Every section of the following code can be copy-pasted into a working GAP session
##  The output of a command is given after a single '#' 
##  A comment is given after a double '#' 

###  
##  In this section, we construct a 24-dimensional representation 
##  for < O8+(2), ρ > over 𝔽_2, where ρ is the triality automorphism,
##  by inducing it from an 8-dimensional representation for O8+(2).

##  p1, p2  are generators of O8+(2) 

p1 := [[1,0,1,1,1,1,0,1],
       [0,0,1,0,1,1,0,1],
       [0,0,1,0,1,0,0,1],
       [0,0,0,0,1,0,0,1],
       [0,0,0,0,1,0,0,0],
       [0,1,0,0,1,0,0,0],
       [0,1,0,0,0,0,1,0],
       [0,0,0,1,0,0,1,0]]*Z(2)^0;;

p2 := [[1,1,0,0,0,0,0,0],
       [0,1,0,0,0,0,0,0],
       [0,0,1,0,0,0,0,0],
       [0,0,0,1,0,0,0,0],
       [0,0,0,1,1,0,0,0],
       [0,0,0,0,1,1,0,0],
       [0,0,0,0,0,0,1,0],
       [1,1,0,0,0,0,0,1]]*Z(2)^0;;

##  s1, s2  are other generators of O8+(2)  ( images of p1, p2 under ρ )

s1 := [[0,0,1,0,0,0,0,0],
       [1,1,0,1,0,0,0,0],
       [1,1,0,0,0,0,0,0],
       [0,1,0,0,0,0,0,0],
       [1,1,0,1,0,0,1,0],
       [1,0,1,0,0,0,0,1],
       [1,0,0,1,1,0,0,1],
       [0,0,0,1,1,1,0,0]]*Z(2)^0;;

s2 := [[1,0,0,0,0,0,0,0],
       [0,1,0,0,0,0,1,0],
       [1,0,1,0,0,1,0,0],
       [0,0,0,1,0,0,0,0],
       [0,1,0,0,1,0,0,0],
       [1,0,0,0,0,1,0,0],
       [0,0,0,0,0,0,1,0],
       [0,0,0,0,0,0,0,1]]*Z(2)^0;;

G := Group(p1,p2);                                                 ##  O8+(2)  ( 8-dimensional representation )

rho8 := GroupHomomorphismByImages( G, G, [p1,p2], [s1,s2] );;      ## the triality automorphism of O8+(2) ( abstract )

## Constructing the induced 24-dimensional representation of O8+(2):3 

rho := KroneckerProduct( PermutationMat( (1,2,3), 3, GF(2) ), 
                         IdentityMat( 8, GF(2) ) );;               ## the triality automorphism of O8+(2) ( 24-dim matrix )
Order(rho); # 3

P1 := IdentityMat( 24, GF(2) );; 
P1{[ 1.. 8]}{[ 1.. 8]} := p1;;
P1{[ 9..16]}{[ 9..16]} := s1^rho8;;
P1{[17..24]}{[17..24]} := s1;;

P2:=IdentityMat( 24, GF(2) );; 
P2{[ 1.. 8]}{[ 1.. 8]} := p2;;
P2{[ 9..16]}{[ 9..16]} := s2^rho8;;
P2{[17..24]}{[17..24]} := s2;;

O8 := Group( P1,  P2 );;                     ## O8+(2) in its 24-dim representation
Size(O8);   #  174182400                    
O8_3 := Group( P1, P2, rho );;               ## O8+(2):3 in its 24-dim representation
Size(O8_3); #  522547200                    

StructureDescription(O8_3);                 
# "O+(8,2) : C3"                             ## The required structure

C_rho := Centralizer( O8, rho );;           
Size( C_rho );   # 12096                     ## = |G2(2)|  => < rho > is indeed the triality automorphism 

Size( Group( rho, rho^P1, rho^P2 ) ); # 522547200  ## = |O8+(2):3| 

## This confirms that α_S(ρ) ⩽ 3, where  S = O8+(2), ρ = triality

Conj_rho := ConjugacyClass( O8_3, rho );;  
Orbs := OrbitsDomain(  C_rho, Conj_rho );;  ## Orbits of C(ρ) on the conjugacy class ρ^S 
Size(Orbs); # 17                            
List(Orbs,Size);
# [ 63, 56, 1512, 1512, 2016, 2016, 63, 2016, 56, 1512, 56, 1512, 63, 378, 1512, 56, 1 ]

Set( Orbs, Orb -> Size( Group( rho, Representative( Orb ) ) ) );
# [ 3, 9, 12, 24, 27, 72, 324, 3456 ]       ##  Sizes of 2-generated subgroups

## Conclusion: Any two conjugates of ρ generate a {2,3}-subgroup of O8+(2):3
###
