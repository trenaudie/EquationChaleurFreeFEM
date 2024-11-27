Je vais reformater ce document en Markdown tout en préservant sa structure et son contenu mathématique.

# Équation de chaleur

## Introduction

Le problème d'équation de chaleur proposé consiste à déterminer la distribution de chaleur d'une pièce, compte tenu des sources de chaleur et de froid présentes dans la pièce, notamment un radiateur, des murs et une fenêtre ouverte.

On donnera une approximation par calcul numérique de la solution stationnaire pour le champ de température T, pour trois ensembles de conditions initiales différentes.

Ensuite nous déterminerons la solution transitoire pour le champ de température T par calcul numérique.

Nous utiliserons pour ceci le cours sur la méthode des éléments finis, et donc une formulation variationnelle du problème, ainsi que le logiciel de calcul numérique FreeFEM++.

## Problème Stationnaire

Le problème stationnaire est le suivant :

$$ k \Delta T = 0 $$

Notons $\Omega$ la pièce et $\Gamma$ les bords de la pièce. $\Gamma = \Gamma_1 + \Gamma_2 + \Gamma_3$ où $\Gamma_1$ est les murs de la pièce, $\Gamma_2$ désigne les bords du radiateur, et $\Gamma_3$ désigne les bords de la fenêtre.

D'où :

$$ \int_\Omega v k \Delta T = 0 $$

D'après la formulation de Green pour le champ de température T, qui est $C^1$:

$$ \int_\Omega v \Delta T = \int_\Omega \nabla v \cdot \nabla T dx + \int_\Gamma v(x) \partial_n T dx $$

On a l'équation finale suivante :

$$ \int_\Omega k \nabla v \cdot \nabla T dx + \int_\Gamma k v(x) \partial_n T dx = 0 $$

## Cas 1 : Flux de chaleur au mur

On considère les conditions suivantes au bord:
- Bord de la fenêtre : Température imposée T = -2°C
- Mur : Flux de chaleur imposée Φ = k∇T·n = −0.31W/m³
- Radiateur : Température imposée T = 50°C

On a donc l'équation :

$$ \int_\Omega k \nabla v \cdot \nabla T dx + \int_\Gamma v(x) \Phi dx = 0 $$

On intègre avec le solveur Conjugate Gradient de FreeFem :

```cpp
solve Laplace(u, v) = 
    int2d(Th)(    // The bilinear part
        dx(u)*dx(v) + dy(u)*dy(v)
    )
    + int1d(Th,c1,c2,a,b,d)( flux / kheat * v )
    //+ int1d(Th, c2)( 0.31 / kheat * v )
    + on(arad, u=tempradiator)
    + on(brad, u=tempradiator)
    + on(crad, u=tempradiator)
    + on(drad, u=tempradiator)
    + on(awind, u=tempwindow)
    + on(bwind, u=tempwindow)
    + on(dwind, u=tempwindow);
```

Resultats :

Avec un flux de chaleur nulle au niveau du mur ( $\Phi = 0 $), on obtient le champ stationnaire suivant 
![alt text](figs/cond1_flux_0.png)

Avec un flux de chaleur non nulle au niveau du mur ( $\Phi = 0.31 W/mK $), on obtient le champ stationnaire suivant 
![alt text](figs/cond1_flux_031.png)


# Cas 2 : Condition de Fourier au niveau du radiateur

Les conditions aux bords sont alors :

- Bord de la fenêtre : Température imposée T = -2°C
- Mur : Température imposée T = -2°C
- Radiateur : Condition de Fourier k∇T·n + h(T - Tf) = 0, avec h = 1W/(m°C) et Tf = 50°C

L'équation devient :

$$ \int_\Omega k \nabla v \cdot \nabla T dx + \int_\Gamma v(x) \frac{h}{k} (T - T_f) dx = 0 $$

On introduit le terme bilinéaire dans le solveur :

```cpp
solve Laplace(u, v) = 
    int2d(Th)(    // The bilinear part
        dx(u)*dx(v) + dy(u)*dy(v)
    )
    + int1d(Th,arad,brad,crad,drad)( v * hfourier / kheat * u)
    - int1d(Th,arad,brad,crad,drad)( v * hfourier / kheat * tempradiator )
    + on(awind, u=tempwindow)
    + on(bwind, u=tempwindow)
    + on(dwind, u=tempwindow)
    + on(a, u=tempwindow)
    + on(b, u=tempwindow)
    + on(c1, u=tempwindow)
    + on(c2, u=tempwindow)
    + on(d, u=tempwindow);
```

Résultats :

![alt text](figs/cond2.png)

# Cas 3 : Condition de Fourier au niveau du mur

Les conditions aux bords sont alors :

- Bord de la fenêtre : Température imposée T = -2°C
- Mur : Condition de Fourier k∇T·n + h(T - Tf) = 0, avec h = 1W/(m°C) et Tf = -2°C
- Radiateur : Température imposée T = 50°C

Résultats :

![alt text](figs/cond3.png)
