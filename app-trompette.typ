#import "@local/modale:1.0.0": *

#show: modale.with(
  title: "APP Trompette",
  subtitle: "Résonances d'un guide d'onde cylindrique et conique : rayonnement, pertes visco-thermiques et déformées modales",
  date: datetime(year: 2026, month: 11, day: 9),
  series: "APP",
  keywords: ("acoustique musicale", "matrice de transfert", "impédance de rayonnement", "pertes visco-thermiques", "analyse modale"),
  abstract: [
    On étudie les résonances d'un instrument à vent assimilé à une chaîne de troncs de cône, depuis le cylindre idéal jusqu'au guide dissipatif et rayonnant. Le formalisme des matrices de transfert conduit à une équation transcendante dont les racines complexes sont les pulsations propres ; les déformées modales s'en déduisent en propageant la condition de sortie jusqu'à l'entrée.
  ],
  cover: true,
  cover-mode: 3,
  toc: true,
)

#symbols(
  ($c$, [Célérité du son], [m/s]),
  ($p$, [Pression acoustique], [Pa]),
  ($U$, [Débit acoustique], [m#super[3]/s]),
  ($Z_c$, [Impédance caractéristique], [Pa·s/m#super[3]]),
  ($Z_r$, [Impédance de rayonnement], [Pa·s/m#super[3]]),
  ($k$, [Nombre d'onde], [rad/m]),
  ($a$, [Rayon du guide], [m]),
  ($Delta L$, [Correction de longueur acoustique], [m]),
  ($Gamma$, [Constante de propagation], [m#super[−1]]),
)

= Cylindre unique
== Sans rayonnement
On prend d'abord le cas d'un cylindre simple, ouvert-ouvert, de section constante et sans considérer le rayonnement sortie. On peut montrer que l'équation vérifiée par la pression acoustique $p(x, t)$ dans le cylindre est

$ (partial^2 p)/(partial x^2) = 1/c^2(partial^2 p)/(partial t^2) $

Dont on peut écrire les solutions dans le domaine fréquentiel sous la forme

- $ P(x, omega) = P^+ (omega) e^(i k x
- ) + P^
(omega) e^(i k x) $ <e1>

On en déduit, via l'équation d'Euler linéarisée, que le débit peut s'écrire

- $ Z_c U(x, omega) = P^+(omega) e^(-i k x)P
^-(omega) e^(i k x) $ <e2>

En combinant ces deux équations écrites au point $x = x_2$, on obtient
- $ cases(P^-(omega) = 1/2 (P_2(omega, x_2)Z
_c U_2(omega, x_2)),
P^+(omega) = 1/2(P_2(omega, x_2) + Z_c U_2(omega, x_2)) ) $

On peut donc réécrire la pression au point $x = x_1$ en substituant ces expressions :

$ P_1(omega, x_1) = P_2(omega, x_2) cos(k(x_2 - x_1)) + j Z_c U_2(omega, x_2) sin(k(x_2 - x_1)) $

En effectuant la même chose pour le débit, on obtient

$ U_1(omega, x_1) = j Z_c^(-1)P_2(omega, x_2) sin(k(x_2 - x_1)) + U_2(omega, x_2) cos(k(x_2 - x_1)) $

On peut donc écrire la relation entre les grandeurs acoustiques aux deux points à l'aide de la matrice de transfert suivante :

$ vec(P_1, U_1) = mat(cos(k(x_2 - x_1)), j Z_c sin(k(x_2 - x_1));
j Z_c^(-1) sin(k(x_2 - x_1)), cos(k(x_2 - x_1))) vec(P_2, U_2) $

D'où on déduit l'expression de l'impédance acoustique au point $x_1$ en fonction de celle en $x_2$ :

$ Z_1 = Z_c (Z_2 cos(k(x_2 - x_1)) + j Z_c sin(k(x_2 - x_1)))/(Z_c cos(k(x_2 - x_1)) + j Z_2 sin(k(x_2 - x_1))) $

En particulier, en prenant $x_1 = 0$ et $x_2 = L$, en ne considérant pas le rayonnement en sortie du cylindre ($Z_2 = 0$) et pas d'embouchure ($Z_1 = 0$), on obtient la relation suivante :

#result[
  $ Z_2 = j Z_c tan(k L)  = 0 => k L = n pi => f_n = n c/(2 L) $
]

Puisqu'on ne considère pas de phénomènes dissipatifs, les solutions modales sont de la forme
$ p(x, t) = sum_n Phi _n (x) q_n (t) \
Phi_n (x) = sin(k_n x) quad q_n (t) = a_n cos(omega_n t) + b_n sin(omega_n t)\
k_n = omega_n / c = (n pi) / L $

Dans le cas d'une embouchure, on aurait un débit nul à l'entrée du cylindre, ce qui correspond à une impédance d'entrée infinie. On chercherait donc les maxima de l'impédance d'entrée pour obtenir les fréquences de résonance du système.

== Avec rayonnement
On se place en basses fréquences : #validity[$k a << 1$] avec $a$ le rayon du cylindre et $k$ le nombre d'onde. On considère les cas où la sortie du cylindre est placée dans un écran infini (rayonne dans un demi-espace libre) ou dans un espace libre (rayonne dans l'espace libre). D'après les transparents de Zegrev, en posant l'impédance de rayonnment de la forme $Z_r = R_r + j X_r$, on a, pour un cylindre bafflé,#note(label: "Source")[Transparents de Zegrev ; Chaigne.]
$ Z_(r, a) = (rho_0 S)/c (((k a)^2)/(2) + j (8 k a)/(3 pi)) $
Et pour un cylindre non bafflé,
$ Z_(r, b) = (rho_0 S)/c (((k a)^2)/(4) + j 0.6133 k a) $

En reprenant l'expression de l'impédance d'entrée du cylindre en fonction de celle de sortie, et en prenant $Z_2 = Z_(r, a\/b)$ cette fois-ci, on obtient :

$ j Z_c tan(k_n L) = -R_r - j X_r $

#hypothesis[
  En réalité, d'après les expressions de l'impédance de rayonnement, on remarque que $R_r << X_r$ lorsque $k a << 1$ (le régime basses fréquences dans lequel on se place), et on peut donc considérer le système conservatif approché :
]

$ Z_c tan(k_n L) = -X_(r, a\/b) $

On conserve donc un nombre d'onde réel et on considère que la présence du rayonnement a avant tout un effet réactif, comme une masse ajoutée au système qui abaisse ses fréquences de résonance. Pour le cas bafflé, on a donc $tan(k_n L) = -(8 k a)\/(3 pi)$. En posant $x_n = k_n L$, on a donc

$ tan(x_n) = -alpha x_n quad text("où") alpha = (8 a)/(3 pi L) $

En considérant $a << L$ (ce qui est le cas pour une trompette), on a $alpha << 1$ et on peut donc considérer que $x_n$ est proche des zéros de la fonction tangente, ce qui revient à écrire

$ x_n approx n pi + epsilon_n quad text("avec") epsilon_n << 1 $

On a alors $tan(x_n) = tan(epsilon_n) approx epsilon_n$ et par ailleurs $-alpha x_n approx -alpha n pi$ au premier ordre. On a donc

$ epsilon_n = - (8 a)/(3 pi L) n pi \
- => k_n L = n pi
(8 a)/(3 pi L) n pi \
=> k_n (L + (8 a)/(3 pi)) = n pi $

On retrouve la notion de correction de longueur acoustique, qui est ici de $8 a\/(3 pi)$ pour un cylindre ouvert dans un demi-espace libre. Pour le cas non bafflé, on obtient de la même manière une correction de longueur acoustique de $0.6133 a$.#note(label: "Ordre de grandeur")[Bafflé : $Delta L approx 0.85 a$. Non bafflé : $Delta L approx 0.61 a$.]

#result(label: "Correction de longueur")[
  $ Delta L_"bafflé" = (8 a)/(3 pi), quad Delta L_"non bafflé" = 0.6133 a $
]

De manière générale, on peut écrire les fréquences de résonance du cylindre avec rayonnement comme

$ f_n = f_n^(0) 1/(1 + (Delta L)/L) approx f_n^(0) (1 - (Delta L)/L) $

où $f_n^(0)$ sont les fréquences de résonance du cylindre ouvert-ouvert sans rayonnement considéré et $Delta L$ la correction de longueur acoustique. On observe que la correction de longueur pour le cas bafflé est plus importante que pour le cas non bafflé, ce qui est cohérent avec le fait que le rayonnement est plus efficace dans le cas bafflé.

L'approximation $R_r << X_r$ nous permet de conserver les mêmes expressions algébriques des solutions modales que pour le cas sans rayonnement, mais avec des fréquences de résonance modifiées numériquement.

#takeaway[
  À basse fréquence, le rayonnement agit comme une masse ajoutée : il allonge acoustiquement le tube de $Delta L$ et abaisse l'ensemble des fréquences de résonance.
]

== Phénomènes dissipatifs
Dans les guides d'onde cylindriques, il existe deux mécanismes principaux de dissipation :
- La viscosité aux parois
- La conduction thermique aux parois

Puisque ces phénomènes apparaissent de manière localisée aux parois, le champ de vitesse dépend de la coordonnée transverse.

#hypothesis[
  La théorie complète existe (Kirchhoff), mais on peut considérer une approche simplifiée (Zwikker & Kosten), qui consiste à garder des équations unidimentionnelles en considérant des grandeurs moyennes sur la surface du guide.
]

D'après les transparents de Zegrev (cf. Chaigne), on peut montrer que cela revient à introduire deux nouvelles grandeurs :#note(label: "Source")[Transparents de Zegrev ; Chaigne.]
- 
$Z_v$ l'impédance linéique en série, modélisant l'effet des frottements visqueux
- 
$Y_t$ l'admittance linéique en parallèle, modélisant l'effet de la conduction thermique
Elles sont telles que
$ (dif P)/(dif x) = -Z_v U quad text("et") quad (dif U)/(dif x) = -Y_t P $

On peut les interpréter de la manière suivante : dans le cas idéal (sans pertes), on a d'après l'équation d'Euler linéarisée que
$ (dif P)/(dif x) = -j omega rho_0 / S U $
Cela revient à considérer une impédance imaginaire pure, c'est-à-dire que les variations locales de pression ne servent qu'à accélérer les particules (c'est une action inertielle seulement). Considèrer $Z_v$ avec une partie réelle (positive) revient à ajouter un terme dissipatif, l'effet des frottements visqueux.

De même, l'équation de conservation de la masse nous donne
$ (dif U)/(dif x) = - j omega S / (rho_0 c^2) P $
Encore une fois, on a une impédance imaginaire pure, et considérer $Y_t$ avec une partie réelle positive revient à prendre en compte la dissipation par conduction thermique.

On peut montrer que ces deux grandeurs sont données par les expressions suivantes :
$
Z_v = (j omega rho) / S (1 - 2/(k_v R) (J_1(k_v R))/(J_0(k_v R)))^(-1) \
- Y_t = j omega chi_s S (1 + (gamma-1) 2/(k_t R) (J_1(k_t R))/(J_0(k_t R)))
$

Avec $S$ la section du guide, $R$ son rayon, $J_0$ et $J_1$ les fonction de Bessel d'ordre 0 et 1, et
- le nombre d'onde de diffusion visqueuse 
$k_v$, où $I_v = mu / (rho_0 c)$ avec $mu$ le coefficient de viscosité dynamique.$ k_v = sqrt((-j omega)/(c I_v)) $
- le nombre d'onde de diffusion thermique 
$k_t$, où $I_t = kappa / (rho_0 c C_p)$ avec $kappa$ le coefficient de conductivité thermique et $C_p$ la chaleur spécifique à pression constante.$ k_t = sqrt((-j omega)/(c I_t)) $
- 
$chi_s$ la compressibilité adiabatique du gaz
- 
$gamma$ le rapport des chaleurs spécifiques à pression et volume constants

On remarque que si $mu = 0$ alors $k_v -> infinity$ et $Z_v = j omega S / rho_0$, de sorte qu'on retrouve l'équation d'Euler linéarisée. Il en va de même pour le cas $kappa = 0$ et l'équation de conservation de la masse.

Ce formalisme introduit une constante de propagation complexe $Gamma$ telle que
$ Gamma = sqrt(Z_v Y_t) = j k_c $
Avec $k_c$ le nombre d'onde complexe. Ainsi, on peut écrire la pression et le débit dans le guide comme
- $ P(x, omega) = P^+ (omega) e^(
- Gamma x) + P^
(omega) e^(Gamma x) \
- U(x, omega) = Z_c^(-1) (P^+ (omega) e^(
- Gamma x)P
- ^
(omega) e^(Gamma x)) $
Où $Z_c$ devient également complexe, et est donnée par
$ Z_c = sqrt(Z_v / Y_t) $

On retrouve donc les mêmes expressions algébriques de matrice de transfert et d'impédance d'entrée que pour le cas conservatif, mais avec des grandeurs complexes (nombre d'onde et impédances).

Enfin, on peut explicitement avoir le lien entre les coefficients $mu$ et $kappa$ et la taille des couches limites visqueuse ($delta_v$) et thermique ($delta_t$) considérées :
$ delta_v = 1/(|k_v|) = sqrt(mu / (rho_0 omega)) \
delta_t = 1/(|k_t|) = sqrt(kappa / (rho_0 C_p omega)) $

= Cône unique
== Cas conservatif
On peut montrer (Chaigne) que la matrice de transfert d'un cône conservatif de rayons $r_1$ et $r_2$ est donnée par
#[
#set text(size: 9.5pt)
$

  vec(P_1, U_1) = mat(r_2/r_1 cos(k l) - (sin(k l))/(k x_1), j (rho c)/(pi r_1 r_2) sin(k l);
  (pi r_1 r_2)/(rho c) ( j sin(k l) (1 + 1/(k^2 x_1 x_2)) + cos(k l)/(j k) (1/x_1 - 1/x_2)), r_1/r_2 cos(k l) + (sin(k l))/(k x_2)) vec(P_2, U_2)

$ <mat_cone>
]

On remarque qu'on a deux impédances acoustiques spécifiques, une pour les ondes aller et une pour les ondes retour, qui sont données par
$ Y_c^(plus.minus) = S/(rho c) (1 plus.minus 1/(j k r)) $
Et on a bien que pour $k r >> 1$, on retrouve le cas des ondes planes.

== Cas dissipatif
Pour prendre en compte les pertes visco-thermiques dans le cône, l'idée est de les considérer pour un rayon moyen $r_m = (r_1 + r_2)\/2$. Ainsi, on peut résoudre de la même manière que pour le cas conservatif, mais en appliquant une correction
- de densité effective
- de compressibilité effective
pour ce rayon moyen.

Cela revient à réécrire @mat_cone en remplaçant
- 
$rho$ par $rho_s$ : $ rho_s = rho (1 + (2 sqrt(-j)) / r_v) $
- 
- $chi_s$ par $chi_t$ : $ chi_t = chi_s (1 + (2 (gamma-1) sqrt(-j)) / r_t) $

Avec $r_v$ et $r_t$ indépendants de $x$ calculés à partir de $r_m$ le rayon moyen du cône :

$ r_v = r_m sqrt(omega / (c l_v)) quad text("et") quad r_t = r_v sqrt(P_r) $
Avec $l_v = mu \/ (rho_0 c)$ et $P_r$ le nombre de Prandtl de l'air.

Les effets visco-thermiques changeant la vitesse des particules, le nombre d'onde devient complexe et est donné par
$
-   k = omega sqrt(rho_s chi_t) = omega / c (1 + (sqrt(-j)) / r_v (1 + (gamma-1) sqrt(P_r)))
$ <k_pertes>

De ce fait, les pulsations sont également complexes et sont données par
- $ omega = k c (1 - (sqrt(-j)) / r_v (1 + (gamma-1)) / sqrt(P_r))) $

Ainsi, une décroissance exponentielle apparaît dans les solutions modales, et on peut montrer que l'amortissement est proportionnel à $1\/r_m$.

#takeaway[
  Les pertes visco-thermiques rendent les pulsations propres complexes : chaque mode décroît exponentiellement, avec un amortissement proportionnel à $1\/r_m$.
]

= Mise en cascade de cônes
On a nécessairement continuité de la pression et du débit aux interfaces entre deux cônes (ou cylindres) consécutifs. Ainsi, on peut remonter à l'impédance d'entrée du système total en multipliant les matrices de transfert de chaque section. On peut donc se ramener à l'expression de l'impédance d'entrée en fonction de celle de sortie (qui sont toutes deux des données du problème), et en déduire les fréquences de résonance du système total. Cette méthode est valable que l'on considère ou non les phénomènes dissipatifs, et que l'on considère ou non le rayonnement à la sortie du système.

= Résolution modale
On considère à présent la concaténation de $M$ troncs de cônes. On commence par déterminer, pour chaque $i in \{1, dots, M\}$,
- le rayon moyen du tronc : 
$r_(m,i) = (r_(1,i) + r_(2,i)) / (2)$
- la pente des parois du tronc : 
- $tan theta_i = (r_(2, i)r
_(1, i)) / L_i$ avec $L_i$ la longueur du tronc
- 
$s_(i, 1)$ et $s_(i, 2)$ les distances de l'apex du cône $i$ à, respectivement, l'extrémité amont du tronc, et l'extrémité aval du tronc. On a donc $s_(i, 2) = s_(i, 1) + L_i$, $r_(i, 1) = s_(i, 1) tan theta_i$ et $r_(i, 2) = s_(i, 2) tan theta_i$.

Ensuite, pour chaque tronc de cône, on construit sa matrice de transfert $M_i$ (@mat_cone) en remplacant $rho$ et $chi_s$ par $rho_s$ et $chi_t$, obtenus à l'aide du rayon moyen du tronc. On calcule également l'impédance de rayonnement en sortie du système, pour le cas bafflé ou non.

On chaîne ensuite ces matrices de transfert, ce qui revient à une matrice de transfert globale $M = M_1 M_2 dots M_M$. Pour la sortie, on a $P_N = Z_r U_N$. On fixe arbitrairement $U_N = 1$ donc $P_N = Z_r$. En écrivant $M$ comme la matrice $(2, 2)$ par blocs $(M_(11), M_(12) ; M_(21), M_(22))$, on obient en entrée
$
  vec(P_0, U_0) = mat(M_(11), M_(12); M_(21), M_(22)) vec(Z_r, 1)
$
Ce qui donne un débit en entrée $U_0 = M_(21) Z_r + M_(22)$.
Dans le cas de la trompette en situation de jeu, l'entrée est fermée par les lèvres du musicien, ce qui impose un débit nul ($U_0 = 0$). Les résonances correspondent donc aux pics de l'impédance d'entrée. On peut donc écrire :

#result(label: "Équation caractéristique")[
  $
    D(omega) = M_(21)(omega) Z_r (omega) + M_(22)(omega) = 0
  $
]

Où $D(omega)$ est la composante de débit à l'entrée du système qui vérifie la condition de sortie. C'est l'équation transcendante du système, dont les solutions sont les pulsations propres du système. Ces solutions sont complexes et on peut les écrire
$
  tilde(omega_n) = omega_(n, d) + j sigma_n
$
On en déduit via la relation de dispersion (@k_pertes) les nombres d'onde propres du système $tilde(k_n)$, complexes également.

Pour obtenir les expressions des déformées modales, il faut réaliser qu'on cherche une solution sans source qui vérifie toutes les conditions aux limites pour $omega = tilde(omega_n)$. Construire $Phi_n$ c'est donc construire un champ de pression qui vérifie quatre conditions :
- l'équation de propagation dans chaque tronc de cône
- la continuité de  $P$ et $U$ aux jonctions
- la condition de sortie  $P_N = Z_r U_N$
- la condition d'entrée  $U_0 = 0$

La condition de sortie fixe le dernier $P_N$, puis on peut propager pour obtenir les $P_i$. Chaque tronc $i$ va de $z_(i-1)$ à $z_i$ (avec $z_0 = 0$ et $z_M = L$). On note aussi $X(z) = (P(z), U(z))^T$. On a donc l'algorithme suivant :

#procedure(
  step("A", [État à la sortie], [
    On a simplement $X_N = (Z_r (tilde(omega_n)), 1)^T$.
  ]),
  step("B", [États aux jonctions], [
    Dans le tronc $i$, la matrice $M_i$ relie l'amont à l'aval : $X_(i-1) = M_i X_i$. On remonte donc :
    $
      X_(M-1) = M_M X_M quad X_(M-2) = M_(M-1) X_(M-1) quad dots quad X_0 = M_1 X_1
    $
  ]),
  step("C", [État à l'intérieur d'un tronc], [
    Pour un point $z in [z_(i-1), z_i]$, on utilise la matrice du morceau compris entre $z$ et $z_i$ :
    $
      X(z) = M_i^text("part") (z -> z_i) X_i
    $
    Où, dans le cas des cônes, $M_i^text("part") (z -> z_i)$ est @mat_cone avec $l = z_i - z$, $s_1 = s(z)$ (distance à l'apex du point $z$), $s_2 = s(z_i)$ (distance à l'apex de la fin du tronc), et les rayons sont $r_1 = s(z) tan theta_i$ et $r_2 = s(z_i) tan theta_i$. Les pertes visco-thermiques sont prises en compte via $k$ et $Z_c$ du tronc $i$ calculés en $tilde(omega_n)$ et donc complexes.
  ]),
  step("D", [La déformée], [
    On a donc $Phi_n(z) = P(z)$ la première composante de $X(z)$. On a bien que
- l'équation de propagation est vérifiée dans chaque tronc 
(puisque par construction, $M$ est la solution locale).
- la continuité aux jonctions 
($X_i$ est repris pour passer d'un tronc à l'autre).
- la condition de sortie est vérifiée puisqu'on a imposé 
$X_N = (Z_r, 1)^T$.
- la condition d'entrée est vérifiée, en propageant jusqu'à l'entrée, la seconde composante de 
$X_0$ (le débit) est nulle.
  ]),
)

En pratique, on peut approximer le résultat en ne gardant que la partie réelle de la déformée dans la solution finale. On peut mesurer l'écart à la solution réelle en prenant le ratio des modules de la partie imaginaire et de la partie réelle.

Les normes $N_n = integral S(z) Re(Phi_n)^2 (z) dif z$ se calculent par intégration numérique une fois les déformées déterminées. Les amplitudes modales sont algébriquement connues.