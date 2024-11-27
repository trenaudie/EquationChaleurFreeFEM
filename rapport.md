## Equation de chaleur 

### Introduction 

Le probleme d'equation de chaleur proposé consiste à déterminer la distribution de chaleur d'une pièce, compte tenu des sources de chaleur et de froid presentes dans la piece, notamment un radiateur, des murs et une fenetre ouverte.

On donnera une approximation par calcul numerique de la solution stationnaire pour le champ de temperature T, pour trois ensembles de conditions initiales differentes. 

Ensuite nous determinerons la solution transitoire pour le champ de temperature T par calcul numerique. 

Nous utiliserons pour ceci le cours sur la methode des elements finis, et donc une formulation variationelle du probleme, ainsi que le logiciel de calcul numerique FreeFEM++.

### Probleme Stationnaire 

Le probleme stationnaire est le suivant :

$$ k \laplacien T = 0  $$ 

Notons $ Omega$ la piece et $ Gamma $  les bords de la piece. $ Gamma =Gamma1 + Gamma2 + Gamma3 $ ou $ Gamma1 $ est les murs de la piece,   $ Gamma2 $ designe les bords du radiateur, et $ Gamma3 $ designe les bords de la fenetre.

D'ou 
$$ \int_\omega v k laplacien T = 0  $$ 
D'apres la formulation de Green pour le champ de temperature T, qui est $C^1$:

$$ \int_\omega v laplacien T = \int_\omega grad v . grad T dx +   \int_\Gamma v(x) ∂_n T dx   $$

On a l'equation finale suivante 
$$  \int_\omega k grad v . grad T dx +   \int_\Gamma k v(x) ∂_n T dx = 0  $$ 

### Cas 1 : Flux de chaleur au mur 

On considere les conditions suivantes au bord: 
Bord de la fenêtre Température imposée T = -2◦C ;
Mur Flux de chaleur imposée Phi = k∇T ·n = −0.31W/m3
Radiateur Température imposée T = 50◦C.

On a donc l'equation 

$$  \int_\omega k grad v . grad T dx +   \int_\Gamma v(x) Phi dx = 0  $$ 


On integre avec le solveur Conjugate Gradient de FreeFem 

Code 
solve Laplace(u, v)
    = int2d(Th)(    // The bilinear part
          dx(u)*dx(v)
        + dy(u)*dy(v)
    )
    +  int1d(Th,c1,c2,a,b,d)( flux / kheat * v )
    //  +  int1d(Th, c2)( 0.31 / kheat * v )
    + on(arad, u=tempradiator)
     + on(brad, u=tempradiator)
     + on(crad, u=tempradiator)
     + on(drad, u=tempradiator)

    + on(awind, u=tempwindow)
     + on(bwind, u=tempwindow)
     + on(dwind, u=tempwindow)
    ;

On obtient


