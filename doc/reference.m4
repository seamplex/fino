changequote([!,!])dnl
% Fino reference sheet

dnl This reference sheet is for [Fino](https://www.seamplex.com/fino) esyscmd([!git describe | sed 's/-/./'!]).

Note that Fino works on top of [wasora](https://www.seamplex.com/wasora), so you should also check the [wasora reference sheet](https://www.seamplex.com/wasora/reference.html) also---not to mention the [wasora RealBook](https://www.seamplex.com/wasora/realbook).
See Fino in action at the [Fino case files](https://www.seamplex.com/fino).

include(reference-toc.md)


# Fino keywords

esyscmd([!../wasora/doc/reference.sh parser kw!])


# Mesh keywords

esyscmd([!../wasora/doc/reference.sh parser kw ../wasora/src/mesh!])


# Special input distributions

TBD.


# Boundary conditions

TBD.


# Result functions

TBD.


# Wasora keywords

esyscmd([!../wasora/doc/reference.sh parser kw ../wasora/src!])


# Fino variables

esyscmd([!../wasora/doc/reference.sh init va!])

