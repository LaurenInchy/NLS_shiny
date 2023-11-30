
# Définition des fonctions d'ajustement 

functionchoices <- c( "Arc tangente" = "arct",
                      "Beta"="betaf", 
                      "Gamma"="gammaf", 
                      "Gamma avec constante"="gammafconst", 
                      "Gompertz"="gompertz", 
                      "Gompertz avec constante" = "gompertzconst",
                      "Gompertz avec freinage" = "gompertzfrein",
                      "Gauss"="gauss",                    
                      "Gauss avec constante"="gaussconst",      
                      "Gauss croissante"="gaussmod", 
                     "Lineaire"="lin", 
                    "Lineaire segmente" = "linseg",
                    "Logaritme 3 parametres"= "log3p", 
                     "Logaritme 4 parametres"="log4p",
                    "Logonormale"="logonorm",
                    "Mitscherlich"="mitscherlich", 
                    "Richards" = "rich",
                    "Richards avec freinage" = "richfrein",
                    "Weibull"="weibull"  )


vitunite <- function(var){
  if (var %in% c("psmg","pfmg","qteaumg")) vu <- "mg" 
  if (var =="teneaupct") vu <- "%" 
  if (var =="qteau") vu <- "mg" 
  if (var =="volmm3") vu <- "mm3" 
  return(vu)
}

tauxmax <- function(modele, a, b, c, d, var){
  if (modele == "mitscherlich"){ v <- a/c ; x <- b ; y <- 0 }# Debouche1979
  if (modele == "gompertz"){ v <- a/(exp(1)*c) ; x <- b ; y <- a/exp(1)}# Debouche1979
  if (modele == "gompertzmod"){ v <- a/(exp(1)*c) ; x <- b ; y <- a/exp(1) +d}
  if (modele == "gaussmod"){ v <- 2*a/(sqrt(2*exp(1))*c) ; x <- b + c/sqrt(2) ; y <- 0.394*a}# Debouche1979
  if (modele == "log3p"){ v <- 0.25*a/c ; x <- b ; y <- 0.5 * a}# Debouche1979
  if (modele == "log4p"){ v <-0.25*a/c ; x <- b  ; y <- 0.5 * a + d}
  if (modele == "weibull"){ v <- exp(-(1-1/c))*a*b*c*(( (c-1)/( b*c) )^(1 - 1/c)); x <- NA ; y <- NA }
  if (modele == "betaf"){ v <- a ; x <- c ; y <- NA}
  if (modele == "arct"){ v <- a / (pi * c) ; x <- b ; y <- 0.5 * a}# Debouche1979
  if (modele == "lin"){ v <- a ; x <- c ; y <- NA}
  if (modele == "linseg"){ v <- NA ; x <- NA ; y <- NA}
  if (modele == "gauss"){v <- NA ; x <- b-c ; y <- a*exp(-.05)}
  if (modele == "gaussconst"){v <- NA ; x <- b-c ; y <- a*exp(-.05)}
  if (modele == "gammaf"){ v <- NA ; x <-  uniroot(function(x){a * exp(-c * x) * (x^(((b - 1) - 1) - 1) * ((b - 1) - 1) * (b - 1)) - a * (exp(-c * x) * c) * (x^((b - 1) - 1) * (b - 1)) - (a * (exp(-c * x) * c) * (x^((b - 1) - 1) * (b - 1)) - a *(exp(-c * x) * c * c) * (x^(b - 1)))},c(1,20), tol=1e-10) ; y <- NA}
  if (modele == "logonorm"){ v <- NA ; x <- uniroot(function(x){-(a * (exp(-((log(x) - b)^2)/(2 * c^2)) * (2 * (1/x * (1/x) -   1/x^2 * (log(x) - b))/(2 * c^2)) - exp(-((log(x) - b)^2)/(2 * c^2)) * (2 * (1/x * (log(x) - b))/(2 * c^2)) * (2 * (1/x *  (log(x) - b))/(2 * c^2)))/(x) - a * (exp(-((log(x) - b)^2)/(2 * c^2)) * (2 * (1/x * (log(x) - b))/(2 * c^2)))/(x)^2 - (a *  (exp(-((log(x) - b)^2)/(2 * c^2)) * (2 * (1/x * (log(x) - b))/(2 * c^2)))/(x)^2 + a * exp(-((log(x) - b)^2)/(2 *   c^2)) * (2 * (x))/((x)^2)^2))}, c(1,20), tol=1e-10) ; y <- NA}
  return(paste( "Vitesse maximale d'accumulation =", round(v,3), vitunite(var), "/j atteint aux coordonnees x = " , round(x,2), "et y = " , round(y,2) , "." ))

}

duree <- function(modele, a, b, c, d){
  
  if (modele == "mitscherlich"){ du <- 2*c }
  if (modele == "gompertz"){ du <- 4*c }
  if (modele == "gompertzfrein"){ du <- c }
  if (modele == "gompertzmod"){ du <- 4*c }
  if (modele == "gaussmod"){ du <- 2 * sqrt(2/pi) * c  }
  if (modele == "log3p"){ du <- 6*c }
  if (modele == "log4p"){ du <- 6*c }
  if (modele == "weibull"){ du <- (log(20)/b)^(1/c) }
  if (modele == "betaf"){ du <- b }
  if (modele == "arct"){ du <- 2*pi*c }
  if (modele == "lin"){ du <- NA }
  if (modele == "linseg"){ du <- NA}
  if (modele == "gauss"){ du <- b }
  if (modele == "gaussconst"){ du <- b }
  if (modele == "logonorm"){ du <- exp(b-c^2)}
  if (modele == "gammaf"){ du <- (b-1)/c}
  
  return(paste( "Temps de la phase de remplissage :", round(du,2) , "jours"))
}

valmax <- function(modele, a, b, c, d, var){
  if (modele == "mitscherlich"){ asym <- a }
  if (modele == "weibull"){ asym <- a }
  if (modele == "gauss"){ asym <- a }
  if (modele == "gaussmod"){ asym <- a }
  if (modele == "logonorm"){asym <- logonorm(a,b,c,x=exp(b-c^2))}
  if (modele == "gammaf"){ y <- k*exp(-(b-1))*(((b-1)/c)^(b-1)) }
  if (modele == "log3p"){ asym <-  a}
  if (modele == "log4p"){ asym <-  a + d}
  if (modele == "arct"){ asym <-  a }
  if (modele == "gompertz"){ asym <-  a }
  if (modele == "gompertzfrein"){ asym <-  a }
  
  return(paste( "Valeur maximale atteinte :", asym , vitunite(var) ))
  
}


#Arctangente d'après Debouche 1979

arct <- function(a, b, c, d=0, x){
  (a / pi) * ( 0.5*pi + atan( (x-b)/c ) )
}


#Beta d'après yin2003

betaf <- function(a, b, c, d=0, x){
  a * x * (  (2*b - c - x)/(2*b - c)  ) * (   ( x / c )^( c / (b - c) )  ) 
}

#gamma d'après Miller 

gammaf <- function(a, b, c, d=0, x){ 
  a * exp(-c*x) * (x^(b-1))
}


#gamma d'après Miller avec constante

gammafconst <- function(a, b, c, d, x){ 
  a * exp(-c*x) * (x^(b-1)) + d
}

#gaussienne d'après Miller

gauss <- function(a, b, c, d=0, x){
  a*exp(-0.5 * ( (x-b)/c )^2) 
}

#gaussienne d'après Miller avec une constante

gaussconst <- function(a, b, c, d, x){
  a*exp(-0.5 * ( (x-b)/c )^2) + d
}

#gaussienne modifiée d'après Deouche

gaussmod <- function(a, b, c, d=0, x){
  a * ( 1 - exp( - ( (x-b) /c )^2 ) ) 
}


#Gompertz d'après Debouche

gompertz <- function(a, b, c, d=0, x){
 a * exp(  -exp( -(x-b) /c )  )
}

#Gompertz avec constante

gompertzconst <- function(a, b, c, d, x){
  d + a * exp(  -exp( -(x-b) /c )  )
}

#Gompertz avec freinage assymetrique

gompertzfrein <- function(a, b, c, d, x){
    a * exp( b*(x-c) - (b/d)*( 1 -  exp( - d *(x-c) ) ) )

}



#L3P d'après Débouche

log3p <- function(a, b, c, d=0, x){
  a / ( 1 + exp( - (x-b) /c ) ) 
}


#L4P mélange Débouche et Robert

log4p <- function(a, b, c, d, x){
  d + ( a /  ( 1 + exp( - (x-b) /c ) )  )
}


#Linéaire

lin <- function(a, b, c=0, d=0, x){
  a*x + b 
}


#Linéaire segmenté : croissance asymptotique

linseg <- function(a,b,c,d,x,y,dta,varform){
  
  # Indiquez a = le nombre de breakpoints
  # b , c , d sont les breakpoints
  
  #pas de forme sur la variance possible
  
  param <- c(b,c,d)  
  eval(parse(text=paste("segmented(lm(", y,"~",x,  ", dta), seg.Z = ~ ", x, ", psi=list(" , x , " = c(", paste(param[1:a], collapse=","), ")))", sep=" " )))
}



#logonormale d'après Miller

logonorm <- function(a , b, c, d=0, x){
  a * exp(-((log(x)-b)^2)/(2*c^2)) / (x+0.000000001)
}



# Mitscherlich d'après Debouche

mitscherlich <- function(a, b, c, d=0, x){
  a * ( 1 - exp( -(x-b) / c ) )
}

# Richards d'après Venus 1978 et Birch 1999

rich <- function(a,b,c,d,x){
  # d doit etre positif
  a*(1+exp(b-c*x))^(-1/d)
}

# Richards avec freinage d'après Welker 1997

richfrein <- function(a, b, c, d, e, x){
  a * ( (b/(b-c)) * exp(-c*d*(x-e)/(b-c)) - (c/(b-c)) * exp(-b*d*(x-e)/(b-c))  )^(-c*(b-c)/d)
}

# Weibull d'après Robert

weibull <- function(a, b, c, d=0, x){
  a * ( 1 - exp (- b * (x^c)) )
}
