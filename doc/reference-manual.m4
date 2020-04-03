## Fino keywords

esyscmd([x../wasora/doc/reference.sh parser kw | sed 's/^#/##/' x])


## Mesh keywords

esyscmd([x../wasora/doc/reference.sh parser kw ../wasora/src/mesh  | sed 's/^#/##/' x])


## Special input distributions

TBD.


## Boundary conditions

TBD.


## Result functions

TBD.


## Wasora keywords

esyscmd([x../wasora/doc/reference.sh parser kw ../wasora/src  | sed 's/^#/##/' x])


## Fino variables

dnl el segundo sed es para reemplazar \text{} por \r{} que es algo de texinfo
esyscmd([x../wasora/doc/reference.sh init va  | sed 's/^#/##/' | sed 's/\\text/\\r/' x])

