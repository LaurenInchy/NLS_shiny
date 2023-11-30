
## Besoins matériels 

library(shiny)
library(shinydashboard)
library(nlme)
library(segmented)

assemblageModel <- function(modele, param){
  paste( modele, "(",  paste(param, collapse=", " ),",", "x)", sep="")
}

## Répertoire des fonctions de croissance
source("FonctionsAjustement.R")


## Données grains
grains <- read.csv("~/Arvalis/Data/dimgr16red_corLauren.csv")
grains$Serie <- as.factor(grains$Serie)
grains$trait <- as.factor(grains$trait)
grains$type <- as.factor(substr(as.character(grains$trait), 2, 2))
grains$grain <- as.factor(grains$grain)
datatemoin <- grains[grains$type=="0",]

colsamp <- c("orange","red", "green", "blue", "pink")


## Application shiny 
ui <- shinyUI(
  
  dashboardPage( skin="yellow", 
                 dashboardHeader(title="Fonctions d'ajustement - remplissage du grain de blé",  titleWidth = "500px"),
                 dashboardSidebar(width = 1),
                 
                 dashboardBody(  
                   
                   #titre première partie
                   fluidRow(
                     box("Visualisation des courbes, estimation des paramètres d'initialisation :", background ="yellow" , width=12)
                   ),
                   
                     fluidRow( 
                       #colonne de gauche - choix modèle, variable, paramètres, zoom
                       column( box(
                             selectInput("modele",
                                         "Modèle :",
                                         choices=functionchoices,
                                         multiple = F),
                             selectInput("variable",
                                         "Variable :",
                                         choices=c("Poids sec (mg)" = "psmg", 
                                                   "Poids frais (mg)" = "pfmg", 
                                                   "Teneur en eau (%)" = "teneaupct",
                                                   "Quantité d'eau (mg)" = "qteaumg", 
                                                   "Volume (mm3)" = "volmm3" ),
                                         multiple = F),
                             width=NULL),
                             
                             fluidRow( 
                               column( 
                                 box(numericInput("a",
                                                  "Paramètre a:",
                                                  value = 75, step=5),
                                     numericInput("b",
                                                  "Paramètre b:",
                                                  value = 1, step=0.01), width=NULL), width=6),
                               column(
                                 box(numericInput("c",
                                                  "Paramètre c:",
                                                  value = 52, step=0.01),
                                     numericInput("d",
                                                  "Paramètre d:",
                                                  value = 0.1, step=0.01),  width=NULL), width=6),  
                               width=4), 
                             box(
                               sliderInput("xlim", 
                                           "Zoom sur X :", 
                                           min=0, max=60, value=c(0,55), step=1),
                               sliderInput("ylim", 
                                           "Zoom sur Y :", 
                                           min=0, max=100, value=c(0,50), step=1),
                               width=NULL
                             ), 
                             width=4), 
                           
                       #colonne de droite : choix des séries, modèle mathématique, visualisation de la courbe
                           column( 
                             box(
                               selectInput("Series",
                                           "Quels lots témoins choisit-on ?",
                                           choices=c("1", "2", "3", "4", "5"),
                                           multiple = T, 
                                           selected=c("1", "2", "3", "4", "5")), 
                               width=NULL),
                             box(
                               solidHeader = T,
                               verbatimTextOutput("modele"),
                               width=NULL),
                             box(
                               solidHeader = T,
                               plotOutput("visucourbe"),
                               width=NULL), 
                             width=8)
                           
                     ),
                   
                   #titre deuxième partie
                     fluidRow(
                       box("Estimation du modèle par NLS ou GNLS à partir des paramètres d'initialisation choisis :", background ="yellow" , width=12)
                     ),
                   
                   
                     fluidRow(
                       
                       #colonne de gauche - modèle
                       column(
                         box(verbatimTextOutput("nlsres"),
                             width=NULL),
                         box("Intervalle de confiance et importance de cette intervalle (taille de l'intervalle pondérée par valeur prédite)",
                             verbatimTextOutput( "confinter"),
                             width=NULL),
                         box("AIC du modèle estimé par Moindres Carrés Non-Linéaires : ",
                             verbatimTextOutput("AICres"),
                             width=NULL),
                         box("Log likelihood : ",
                             verbatimTextOutput("loglike"),
                             width=NULL),
                         width=6),
                       
                       #colonne de droite - résidus
                       column(
                         box(selectInput("varianceFunction",
                                         "Forme de la variance des résidus :",
                                         choices=c("Homoscedasticite" = "homo",
                                                   "Puissance de Y" = "varPower",
                                                   "Exponentiel de Y" = "varExp",
                                                   "Constante plus puissance de Y" ="varConstPower",
                                                   "Poids fixes en fonction de Y"= "varFixed")),
                          width=NULL),
                         box( plotOutput("residus"), width=NULL) ,
                         box( "Caractéristiques biologiques du modèle :",
                              verbatimTextOutput("tauxmax"), 
                              verbatimTextOutput("duree"), width=NULL
                         ),
                         width=6)

                     )

  )
  )
)
  
  
  


# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  
  # formule mathématique
  output$modele <- renderPrint({ 
    eval(parse(text=input$modele))   })
  
  # choix de la série
  dta <- reactive({
    datatemoin[datatemoin$Serie %in% input$Series , ]
  })
  
  # visualisation de la courbe sur les données
  output$visucourbe <- renderPlot({ 
    
    x <- seq(0.1,54,2)
    param <- c(input$a, input$b, input$c, input$d)
    
    with(dta(), 
         plot(daa, eval(parse(text=input$variable)), 
              type="p", pch=16, col=colsamp[as.factor(Serie)],
              xlim=input$xlim, ylim=input$ylim,
              xlab = "Jours après anthèse",
              ylab=paste(input$variable)))

     if (input$modele == "linseg"){
       
       nb <- round(param[1])
       param2 <- param
       param2[1] <- 0                             
       param2[nb+2] <- 54 
       dta <- dta()

      for ( i in 1:(nb+1) ){
        segmentdta <- dta[ which( (dta$daa >= param2[i]) &(dta$daa <= param2[i+1] )) , ]
        modlm <- lm( eval(parse(text=input$variable)) ~daa,   segmentdta)$coefficients
        lines( seq(param2[i], param2[i+1], 1) ,  modlm[1]+modlm[2]*seq(param2[i], param2[i+1], 1),
               type="l", lwd=2, col="red")
        if(!is.null(nlsmod())){ plot(nlsmod() , add=T , col="grey", lwd=2) }
      }
       
    
       
       }else{
    
         lines(x, 
               eval(parse(text=assemblageModel(input$modele, param))),
               type="l", lwd=2, col="red") 
         derivexp=parse(text=paste(body(eval(parse(text=input$modele)))[2]))
         firstderiv <- eval(parse(text=paste( "function(a,b,c,d,x){", paste(deparse(D(derivexp, 'x')), collapse="") , "}" )))

         lines(x, eval(parse(text=paste("10*firstderiv(", paste(param, collapse=","), ", x=x)"))), col="blue", lwd=2)
         legend("bottomright", leg=c("Courbe","10 * Dérivée"), col=c("red","blue"), lwd=2 )
         
       }
    
    legend("topleft", input$Series, col=colsamp, pch=16 )
    
    })
  
  nlsmod <- reactive ({
    
    dtaclean <- subset( dta(), !is.na(  eval(parse(text=input$variable )) ) )
    
    aEl <- input$a
    bEl <- input$b
    cEl <- input$c
    dEl <- input$d
    
    if (input$modele=="linseg"){
      
      linseg(aEl, bEl, cEl, dEl, "daa" , input$variable , dtaclean) 
      
    }else{
    
    if (input$varianceFunction == "homo") {

          nls( 
            
            if (formals(input$modele)$d == 0){
              
              if (formals(input$modele)$c == 0){
                
            eval(parse(text=paste(input$variable, "~", input$modele,"(a, b, x=daa)", sep="")))
                }else{
                  eval(parse(text=paste(input$variable, "~", input$modele,"(a, b, c, x=daa)", sep="")))
                  }
            }else{
              eval(parse(text=paste(input$variable, "~", input$modele,"(a, b, c, d, x=daa)", sep="")))
            },
               data=dtaclean,
               start=
                 if (formals(input$modele)$d == 0){
                   
                   if (formals(input$modele)$c == 0){
                     list( a=aEl, b=bEl )
                   }else{
                     list( a=aEl, b=bEl, c=cEl )
                   }
                 }else{
                   list( a=aEl, b=bEl, c=cEl, d=dEl )
                  } ,
            nls.control(maxiter = 500)
    )
    
      
      }else{
        
        modform <- if (formals(input$modele)$d == 0){
          if (formals(input$modele)$c == 0){
            eval(parse(text=paste(input$variable, "~", input$modele,"(a, b, x=daa)", sep="")))
          }else{
            eval(parse(text=paste(input$variable, "~", input$modele,"(a, b, c, x=daa)", sep="")))
          }
        }else{
          eval(parse(text=paste(input$variable, "~", input$modele,"(a, b, c, d, x=daa)", sep="")))
        }
      
        gnls(  modform  ,
          
          
          data=dtaclean,
          
          
          start=
            if (formals(input$modele)$d == 0){
              if (formals(input$modele)$c == 0){
                list( a=aEl, b=bEl )
              }else{
                list( a=aEl, b=bEl, c=cEl )
              }
            }else{
              list( a=aEl, b=bEl, c=cEl, d=dEl )
            } ,
          
          
          weights= 
            if (input$varianceFunction == "varFixed") {
              eval(parse(text=paste(input$varianceFunction, "( ~", input$variable, ")", sep="") ))
            }else{ 
              eval(parse(text=paste(input$varianceFunction, "( form = ~", input$variable, ")", sep="") ))
              },
          
          
          control=gnlsControl(nlsTol=1)
        )
      }
    }
    
  })
  
  output$nlsres <- renderPrint({ 
    summary( nlsmod()   , correlation = T )
  })
  
  aest <- reactive( coefficients( nlsmod() )["a"] )
  best <- reactive( coefficients( nlsmod() )["b"] )
  cest <- reactive( coefficients( nlsmod() )["c"] )
  dest <- reactive( ifelse(formals(input$modele)$d != 0, coefficients( nlsmod() )["d"], 0 ))
  
  
  output$confinter <- renderPrint({
    tempconf <- as.data.frame(confint( nlsmod() ) )
    if( dest() != 0 ){
      tempconf$perc <- abs( 100* abs( tempconf[,2]-tempconf[,1] ) / c( aest(), best(), cest(), dest() ) )
    }else{
      tempconf$perc <- abs( 100* abs( tempconf[,2]-tempconf[,1] ) / c( aest(), best(), cest() ) )
    }
   tempconf
  })
  
  output$AICres <- renderPrint({
  AIC( nlsmod() ) 
  })
  
  output$loglike <- renderPrint({
    logLik( nlsmod() ) 
  })
  
  output$residus <- renderPlot({
    plot( resid( nlsmod() , type="pearson") ~ fitted( nlsmod() ),
         ylab="residus de Pearson", xlab="valeurs prédites", pch=16)
    lines(smooth.spline( fitted( nlsmod()), resid( nlsmod() , type = "pearson")  , df=5 ), lwd=2)
    legend("bottomleft", lwd=2, col="black", legend=c("spline"), bty="n")
    abline(h=0, col="red")

  })
  
 
  
   output$tauxmax <- renderPrint({ 
    print(paste(tauxmax(input$modele, aest(), best(), cest(), dest() , input$variable ), "et ",
          valmax(input$modele, aest(), best(), cest(), dest() , input$variable )))
  })
  
   output$duree <- renderPrint({ 
     print(duree(input$modele, aest(), best(), cest(), dest() ))
   })
   
})

# Run the application 
shinyApp(ui = ui, server = server)

