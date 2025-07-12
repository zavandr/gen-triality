###
##  This code for GAP ( version 4.14.0 of 2024-12-05 ) accompanies the paper
##  
##  "On generation by triality automorphisms"
## 
##  by Danila O. Revin and Andrei V. Zavarnitsine
##
##  Date: 11 July, 2025

##  Every section of the following code can be copy-pasted into a working GAP session.
##  The output of a command is given after a single '#'. 
##  A comment is given after a double '#'. 

###  
##  Section 1.  S = O8+(2)
##
##  In this section, we construct a 24-dimensional representation 
##  for < S, ρ > over 𝔽_2, where ρ is a triality automorphism of S,
##  by inducing it from an 8-dimensional representation for S.
##  Then we prove that α_S(ρ) = 3 by finding the orders of subgroups generated 
##  by a pair of elements conjugate to ρ. We also show that α_S(τ) = 2, 
##  where τ ∊ Sρ has order 3 and is not conjugate to ρ.

##  x and y  are generators of O8+(2) : 

x := [[1,0,1,1,1,1,0,1],
      [0,0,1,0,1,1,0,1],
      [0,0,1,0,1,0,0,1],
      [0,0,0,0,1,0,0,1],
      [0,0,0,0,1,0,0,0],
      [0,1,0,0,1,0,0,0],
      [0,1,0,0,0,0,1,0],
      [0,0,0,1,0,0,1,0]]*Z(2)^0;;

y := [[1,1,0,0,0,0,0,0],
      [0,1,0,0,0,0,0,0],
      [0,0,1,0,0,0,0,0],
      [0,0,0,1,0,0,0,0],
      [0,0,0,1,1,0,0,0],
      [0,0,0,0,1,1,0,0],
      [0,0,0,0,0,0,1,0],
      [1,1,0,0,0,0,0,1]]*Z(2)^0;;

##  xr and yr  are other generators of O8+(2)  ( the images of x and y under ρ )

xr := [[0,0,1,0,0,0,0,0],
       [1,1,0,1,0,0,0,0],
       [1,1,0,0,0,0,0,0],
       [0,1,0,0,0,0,0,0],
       [1,1,0,1,0,0,1,0],
       [1,0,1,0,0,0,0,1],
       [1,0,0,1,1,0,0,1],
       [0,0,0,1,1,1,0,0]]*Z(2)^0;;

yr := [[1,0,0,0,0,0,0,0],
       [0,1,0,0,0,0,1,0],
       [1,0,1,0,0,1,0,0],
       [0,0,0,1,0,0,0,0],
       [0,1,0,0,1,0,0,0],
       [1,0,0,0,0,1,0,0],
       [0,0,0,0,0,0,1,0],
       [0,0,0,0,0,0,0,1]]*Z(2)^0;;

G := Group(x,y);;                            ##  O8+(2)  ( 8-dimensional representation )

StructureDescription(G); # "O+(8,2)"         ##  < x, y > is indeed O8+(2)

xr in G; # true                              ##  confirming that
yr in G; # true                              ##  xr, yr ∊ < x, y >  

rho := GroupHomomorphismByImages( 
       G, G, [ x, y ], [ xr, yr ] );;        ## a triality automorphism of O8+(2) ( abstract )

rho <> fail;  #  true                        ## confirming that the map  [ x, y ] -> [ xr, yr ] 
                                             ## indeed extends to an automorphism of O8+(2)                                            

Order(rho); #  3                             ## |ρ| = 3   ( as an abstract automorphism )

## Constructing the induced 24-dimensional representation of O8+(2):3 :

rho0 := KroneckerProduct(                    ## triality automorphism ρ of O8+(2) 
  PermutationMat( (1,2,3), 3, GF(2) ),       ## ( a 24-dimensional permutation matrix ) 
  IdentityMat( 8, GF(2) ) );;                

Order(rho0); # 3                             ## |ρ| = 3   ( as a 24-dimensional matrix )

x0 := IdentityMat( 24, GF(2) );; 
x0{[ 1.. 8]}{[ 1.. 8]} := x;;
x0{[ 9..16]}{[ 9..16]} := xr^rho;;
x0{[17..24]}{[17..24]} := xr;;

y0:=IdentityMat( 24, GF(2) );; 
y0{[ 1.. 8]}{[ 1.. 8]} := y;;
y0{[ 9..16]}{[ 9..16]} := yr^rho;;
y0{[17..24]}{[17..24]} := yr;;

O8 := Group( x0,  y0 );;                     ## O8+(2) in its 24-dimensional representation
Size(O8);   #  174182400                     ## = |O8+(2)|

O8_3 := Group( x0, y0, rho0 );;              ## O8+(2):3 in its 24-dimensional representation
Size(O8_3); #  522547200                     ## = |O8+(2):3|

StructureDescription(O8_3); 
# "O+(8,2) : C3"                             ## the required structure

C_rho := Centralizer( O8, rho0 );;           
Size( C_rho );   # 12096                     ## = |G2(2)|  =>  ρ is indeed a triality automorphism 

Group( rho0, rho0^x0, rho0^y0 ) = O8_3; 
# true                                       ## This confirms that α_S(ρ) ⩽ 3, where  S = O8+(2)

Conj_rho := ConjugacyClass( O8_3, rho0 );;   ## the conjugacy class ρ^S 

Orbs := OrbitsDomain(  C_rho, Conj_rho );;   ## orbits of C(ρ) on ρ^S
Size(Orbs); # 17                             ## number of orbits

List(Orbs,Size);                             ## sizes of orbits
# [ 63, 56, 1512, 1512, 2016, 2016, 63, 2016, 
#  56, 1512, 56, 1512, 63, 378, 1512, 56, 1 ]

Set( Orbs, Orb -> Size( 
   Group( rho0, Representative( Orb ) ) ) );
# [ 3, 9, 12, 24, 27, 72, 324, 3456 ]        ## sizes of 2-generated subgroups  < ρ, ρ^g >, g ∊ S

## Conclusion: Any two conjugates of ρ generate a {2,3}-subgroup of O8+(2):3
##             In particular, α_S(ρ) = 3, where  S = O8+(2)

##  We now prove that α_S(τ) = 2, where τ is a nontriality, 
##  i.e. τ ∊ Sρ has order 3 and τ is not conjugate to ρ

tau := ( rho0 * x0^3 )^4;;                   ## τ 

Order( tau );  # 3                           ## |τ| = 3

Size( Centralizer( O8, tau ) ); # 216        ## = |PGU(3,2)|  =>  τ is the required nontriality automorphism

Size( Group( tau, tau^x0 ) ); # 522547200    ## = |O8+(2):3|  =>  α_S(τ) ⩽ 2
                                             ## Clearly, we also have  α_S(τ) ⩾ 2 

## Conclusion: α_S(τ) = 2, where S = O8+(2), τ = nontriality 
###

###
## Section 2. S = O8+(3)

##  We prove that α_S(ρ) = 3, where ρ is a triality automorphism of S 
##  by finding the orders of subgroups generated by a pair of elements conjugate to ρ. 
##  We also show that α_S(τ) = 2, where τ ∊ Sρ has order 3 and is not conjugate to ρ.

O8_S4 := AtlasGroup( "O8+(3).S4" );;         ## extension O8+(3):S_4

Gens := GeneratorsOfGroup( O8_S4 );;  
Size( Gens ); # 2                            ## there are two generators of O8+(3):S_4

x0 := Gens[1];; y0 := Gens[2];;              ## generators
Order(x0); # 24
Order(y0); # 20                        

O8_A4 := DerivedSubgroup(O8_S4);; 
O8_K4 := DerivedSubgroup(O8_A4);;
O8    := DerivedSubgroup(O8_K4);;            ## O8+(3) as the third derived subgroup of O8+(3):S_4
Size(O8); # 4952179814400                    ## = |O8+(3)|

## Finding generators of O8+(3) :

x := x0^4;;   
y := y0^4;;  

x in O8; # true
y in O8; # true

Size( Group( x, y ) ); # 4952179814400       ## = |O8+(3)|  => x and y generate O8+(3)
                                    
rho := (x0*y0)^4 ;;                          ## ρ
Order( rho ); # 3                            ## |ρ| = 3
rho in O8;  # false                          ## ρ is outside O8+(3)

C_rho := Centralizer( O8, rho );;            ## C_S(ρ), where S = O8+(3)
Size( C_rho );  #  4245696                   ## = |G2(3)|  =>  ρ is indeed a triality automorphism of O8+(3)

O8_3 := Group( rho, rho^x, rho^y );;
Size( O8_3 ); #  14856539443200              ## = |O8+(3):3|  =>  this is indeed O8+(3):3
                                             ## This also confirms that  α_S(ρ) ⩽ 3, where S = O8+(3)

Conj_rho := ConjugacyClass( O8_3, rho );;
Orbs := OrbitsDomain(  C_rho, Conj_rho );;   ## orbits of C(ρ) on the conjugacy class ρ^S 
Size( Orbs ); # 18
List( Orbs, Size );                          ## sizes of orbits
# [ 78624, 176904, 157248, 78624, 78624, 
#  176904, 176904, 176904, 44226, 17472, 
#  728, 728, 351, 728, 351, 728, 351, 1 ]

Set( Orbs, Orb -> Size( 
     Group( rho, Representative( Orb ) ) ) );
# [ 3, 9, 12, 24, 27, 72, 243, 324, 3456 ]   ## sizes of 2-generated subgroups < ρ, ρ^g >, g ∊ S

## Conclusion: Any two conjugates of ρ generate a {2,3}-subgroup of O8+(3):3
##             In particular, α_S(ρ) = 3, where  S = O8+(3)


##  We now prove that α_S(τ) = 2, where τ is a nontriality, 
##  i.e. τ ∊ Sρ has order 3 and τ is not conjugate to ρ.

tau := (rho * x^4 * y)^4;;                   ## τ

Order(tau);   # 3                            ## |τ| = 3

C_tau := Centralizer(O8, tau);;
Size(C_tau); # 5832                          ## = |[3^5].SL(2,3)|  =>  τ is the required nontriality automorphism

Size( Group( tau, tau^x ) ); 
# 14856539443200                             ## = |O8+(3):3|  =>  α_S(τ) ⩽ 2
                                             ## Clearly, we also have  α_S(τ) ⩾ 2 

## Conclusion: α_S(τ) = 2,  where S = O8+(3)
###

### END ###
