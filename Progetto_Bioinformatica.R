### PROGETTO BIOINFORMATICA A.A.2019-2020 
## STUDENTI: 
# BELLEMO FRANCESCA IN1000277
# GATTAFONI DANIELE IN1000268
# MINISINI VALENTINA IN1000250
# SACCOMANO GIULIA IN1000251

setwd("/Users/giuli_pet/Desktop/ProgettoBioinformatica/")
library(RMySQL)
library(dplyr)
library(tidyr)

# Definizione dei parametri di input:

# parametro di testo che viene riportato in ogni riga dell'output, sotto la colonna 'Gruppo'
gruppo = c("F.Bellemo, D.Gattafoni, V.Minisini, G.Saccomano") 
# file di testo contenente un elenco di identificativi da ricercare
file_id = "input.txt"
# parametro che permetta di limitare la ricerca ad uno solo dei valori tra i 53 disponibili nel campo ExpScores della tabella gtextransExpr.csv
number = 25
# nome del file di output dove saranno scritti i risultati della funzione
file_out = "Risultato.txt"

# Creazione funzione commentata:
funzione_ProgettoBioinformatica <- function(gruppo = "Gruppo", file_id, input_number = 1, file_out = "Risultato.txt"){
  
  ### PARTE A
  # Creazione del nome del file di testo contenente i Log
  name_log <- paste0("Progetto_Bioinformatica.",format(Sys.time(),"%d-%m-%Y.%H-%M-%S"),".Log.txt")
  
  # variabile di scrittura del file di Log
  out_log <- file(name_log,"w")
  
  # Il file di Log contiene:
  # elenco degli input forniti
  write("ELENCO INPUT FORNITI:",out_log)
  
  # Esecuzione controllo sul tipo e la lunghezza dei valori di input:
  if(is.character(gruppo) && length(gruppo)==1 &&
     is.character(file_id) && length(file_id)==1 &&
     is.numeric(input_number) && length(input_number)==1 && input_number>0 && input_number<54 && 
     is.character(file_out) && grepl(".+\\.txt", file_out) && length(file_out)==1 ){
    
    write(paste0("gruppo=",gruppo),out_log,append=TRUE)
    write(paste0("file_id=",file_id),out_log,append=TRUE)
    write(paste0("input_number=",input_number),out_log,append=TRUE)
    write(paste0("file_out=",file_out),out_log,append=TRUE)
    
    # l'elenco degli oggetti prodotti dalla funzione
    write("\nELENCO OGGETTI PRODOTTI:",out_log,append=TRUE)
    
    # Controllo esistenza del "file_id"
    if(length(file_id)==1 && grepl(".+\\.txt$",file_id)){
      
      # Se file non esistente, avviene la rilevazione dell'errore e la segnalazione nel file di Log
      if(!file.exists(file_id)){
        write("\nERRORI:",out_log,append=TRUE)
        write(paste(file_id,"non trovato"),out_log,append=TRUE)
        close(out_log)
        stop(paste(file_id,"non trovato"))
      }
      
      # lettura ed estrazione degli identificativi nel file di input
      read_id <- read.table(file_id, header = FALSE) 
      input_id <- as.character(read_id[,"V1"]) 
      
      # apertura connessione su UCSC Genome Browser nell'assembly di Uomo HG38
      mydb <- dbConnect(MySQL(), user="genome", host="genome-mysql.soe.ucsc.edu", port=3306, dbname="hg38")
      
      # importazione della tabella refSeqFuncElems tramite file esterno, utilizzata nella Parte C
      refSeqFuncElems <- read.csv("refSeqFuncElems.csv",header=TRUE,sep=";")
      # scrittura dell'oggetto nel file di Log
      write("refSeqFuncElems",out_log,append=TRUE) 
      
      # importazione della tabella gtexTranscExpr tramite file esterno, utilizzata nella Parte D
      tabella_gtex <- read.csv("gtextransExpr.csv", header = TRUE, sep = ';')
      # scrittura dell'oggetto nel file di Log
      write("tabella_gtex",out_log,append=TRUE)
      
      # estrazione della colonna expScores  
      all_expScores <- data.frame(tabella_gtex[,"expScores"])
      # costruzione di una tabella di soli expScores suddivisi singolarmente
      expScores <- separate(all_expScores, colnames(all_expScores), as.character(c(1:53)), sep = ",")
      
      # inizializzazione dei files tibble relativi alle parti B, C, D
      final_tibble_B <- tibble()
      final_tibble_C <- tibble()
      final_tibble_D <- tibble()
      
      # per ogni identificativo nel file input_id:
      for (identificativo in input_id){
        
        ### PARTE B
        # estrazione valori dei campi geneId, geneName, transcriptId ottenuti dalla ricerca rispetto:
        # Primary Table wgEncodeGencodeAttrsV27, Database hg38, Field geneName e/o geneId
        prima_ricerca <- dbGetQuery(mydb,paste0("SELECT geneId, geneName, transcriptId FROM wgEncodeGencodeAttrsV27	
                                          WHERE geneId ='",identificativo,"' OR geneName = '",identificativo,"';"))
        
        # calcolo numero identicativi corrispondenti alla prima ricerca
        num <- nrow(prima_ricerca)
        
        # in caso di zero corrispondenze:
        if (num == 0){
          # costruzione righe come data.frame per i tibble finali popolando i campi con NA
          colonne_risultato_B <- data.frame(gruppo = gruppo, identificativo = identificativo, geneId = NA, geneName = NA, transcriptId = NA, 
                                            rnaAcc = NA, chrom = NA, strand = NA, txStart = NA, txEnd = NA)
          colonne_risultato_C <- data.frame(identificativo = NA, name = NA, name = NA, soTerm = NA, experiment = NA, function_ = NA)
          colonne_risultato_D <- data.frame(gruppo = gruppo, identificativo = identificativo, name = NA, name2 = NA, expScore = NA, rapporto = NA)
          
          # aggiunta righe ai tibble finali
          final_tibble_B <- bind_rows(final_tibble_B, colonne_risultato_B)
          final_tibble_C <- bind_rows(final_tibble_C, colonne_risultato_C)
          final_tibble_D <- bind_rows(final_tibble_D, colonne_risultato_D)
          
        } else {
          # per ogni corrispondenza della prima ricerca: 
          for (k in 1:num){
            
            # estrazione dei valore dei campi geneId, geneName e transcriptId 
            geneId <- prima_ricerca[k, "geneId"]
            geneName <- prima_ricerca[k, "geneName"]
            transcriptId <- prima_ricerca[k, "transcriptId"]
            
            # estrazione degli eventuali valori dei campi rnaAcc ottenuti dalla ricerca degli identificativi transcriptId rispetto:
            # Primary Table wgEncodeGencodeRefSeqV27, Database hg38, Field transcriptId
            seconda_ricerca <- dbGetQuery(mydb, paste0("SELECT rnaAcc FROM wgEncodeGencodeRefSeqV27 
                                                       WHERE transcriptId = '",transcriptId,"';"))
            
            # in caso di zero corrispondenze della seconda_ricerca:
            if (nrow(seconda_ricerca) == 0){
              # popolazione di tutti i campi della prima_ricerca e NA per la seconda_ricerca
              colonne_ricerca2 <- data.frame(gruppo = gruppo, identificativo = identificativo, geneId = geneId, 
                                             geneName = geneName, transcriptId = transcriptId, rnaAcc = NA)
            } else {
              # in caso di corrispondenza della seconda_ricerca, estrazione del valore del campo rnaAcc
              rnaAcc <- seconda_ricerca[, "rnaAcc"]
              # popolazione di tutti i campi della prima_ricerca e della seconda_ricerca
              colonne_ricerca2 <- data.frame(gruppo = gruppo, identificativo = identificativo, geneId = geneId, 
                                             geneName = geneName, transcriptId = transcriptId, rnaAcc = rnaAcc)
            }
            
            # estrazione della posizione genomica (chrom, strand, txStart, txEnd) degli identificativi transcriptId rispetto:
            # Primary Table wgEncodeGencodeBasicV27, Database hg38, Field name
            terza_ricerca <- dbGetQuery(mydb, paste0("SELECT chrom, strand, txStart, txEnd FROM wgEncodeGencodeBasicV27 
                                                     WHERE name = '",transcriptId,"';"))
            
            # in caso di zero corrispondenze della terza_ricerca:
            if (nrow(terza_ricerca) == 0){
              # popolazione di tutti i campi delle ricerche successive con NA
              colonne_ricerca3 <- data.frame(chrom = NA, strand = NA, txStart = NA, txEnd = NA)
              colonne_ricerca4 <- data.frame(name = NA)
              colonne_ricerca5 <- data.frame(name = NA, soTerm = NA, experiment = NA, function_ = NA)
            } else {
              # in caso di corrispondenza della terza_ricerca, estrazione dei valori di chrom, strand, txStart, txEnd
              chrom <- terza_ricerca[, "chrom"]
              strand <- terza_ricerca[, "strand"]
              txStart <- terza_ricerca[, "txStart"]
              txEnd <- terza_ricerca[, "txEnd"]
              
              # popolazione di tutti i campi della terza_ricerca
              colonne_ricerca3 <- data.frame(chrom = chrom, strand = strand, txStart = txStart, txEnd = txEnd)
              
              ### PARTE C
              # estrazione dei valori del campo name ottenuti selezionando i dati compresi nelle posizioni genomiche 
              # relative alla terza_ricerca rispetto:
              # Primary Table ucscGenePfam, Database hg38
              quarta_ricerca <- dbGetQuery(mydb, paste0("SELECT name FROM ucscGenePfam WHERE chrom = '",chrom,"' AND strand = '",strand,"' 
                                          AND chromStart >= '",txStart,"' AND chromEnd <= '",txEnd,"';"))
              
              # in caso di zero corrispondenze della quarta_ricerca:
              if (nrow(quarta_ricerca)==0) {
                # popolazione del campo name con NA
                colonne_ricerca4 <- data.frame(name = NA)
              } else {
                # in caso di corrispondenza della quarta_ricerca, estrazione del valore name
                name <- quarta_ricerca[, "name"]
                # popolazione del campo name della quarta_ricerca
                colonne_ricerca4 <- data.frame(name = name)
              }
              
              # ricerca degli indici di riga degli eventuali valori dei campi name, soTerm, experiment, function ottenuti selezionando i dati 
              # compresi nelle posizioni genomiche rispetto:
              # Tabella refSeqFuncElems (tabella costruita in formato csv a partire dall'annotation track RefSeqFuncElems
              # Big Bed File: /gbdb/hg38/ncbiRefSeq/refSeqFuncElems.bb)
              indici_estrazione <- subset(refSeqFuncElems, X.chrom == chrom & strand == strand & thickStart >= txStart & thickEnd <= txEnd)
              
              # in caso di zero corrispondenze della ricerca indici_estrazione:
              if (nrow(indici_estrazione)==0) {
                # popolazione di tutti i campi della ricerca con NA
                colonne_ricerca5 <- data.frame(name = NA, soTerm = NA, experiment = NA, function_ = NA)
              } else {
                # in caso di corrispondenza della ricerca, estrazione dei valori name, soTerm, experiment, function.
                name <- refSeqFuncElems[indici_estrazione,"name"]
                soTerm <- refSeqFuncElems[indici_estrazione,"soTerm"]
                experiment <- refSeqFuncElems[indici_estrazione,"experiment"]
                function. <- refSeqFuncElems[indici_estrazione,"function."]
                
                # popolazione del data.frame con i campi name, soTerm, experiment, function. della ricerca
                colonne_ricerca5 <- data.frame(name = name,soTerm = soTerm,experiment = experiment,function_ = function.)
              }
            }
            # costruzione righe come data.frame per i tibble finali B,C
            colonne_risultato_B <- data.frame(colonne_ricerca2, colonne_ricerca3)
            colonne_risultato_C <- data.frame(identificativo = transcriptId, colonne_ricerca4, colonne_ricerca5)
            
            # aggiunta righe ai tibble finali B,C
            final_tibble_B <- bind_rows(final_tibble_B, colonne_risultato_B)
            final_tibble_C <- bind_rows(final_tibble_C, colonne_risultato_C)
            
            # troncamento della stringa transcriptId prima del punto (.) per ottenere corrispondenze nella ricerca per la Parte D
            transcriptId_D <- substr(transcriptId, 0, 15)
            
            ### PARTE D
            # ricerca degli indici di riga degli eventuali sinonimi dell'identificativo transcriptId_D ottenuti rispetto:
            # Tabella gtexTranscExpr.csv
            indice_riga <- grep(transcriptId_D, tabella_gtex[,"name"], perl = TRUE, value = FALSE)
            
            # in caso di zero corrispondenze della ricerca indice_riga:
            if(length(indice_riga)==0){ 
              # popolazione di tutti i campi della corrispondente e della successiva ricerca con NA
              colonne_ricerca_D <- data.frame(name = NA, name2 = NA)
              colonna_ricerca_D2 <- data.frame(expScore = NA, rapporto = NA)
            } else {
              # estrazione dei valori name, name2 
              name1 <- tabella_gtex[indice_riga,"name"]
              name2 <- tabella_gtex[indice_riga,"name2"]
              
              # popolazione del data.frame con i campi name, name2 della ricerca
              colonne_ricerca_D <- data.frame(name = name1, name2 = name2)
              
              # estrazione dell'eventuale valore di score del campo selezionato di expScores per l'identificativo
              column_expScores <- as.numeric(expScores[,input_number])
              expScore <- column_expScores[indice_riga]
              
              # calcolo valor medio di tutti gli score dello stesso campo
              valore_medio <- mean(column_expScores)
              # calcolo rapporto del valore di expScores rispetto al valor medio di tutti gli score dello stesso campo
              rapporto <- expScore/valore_medio
            
              # popolazione del data.frame con i campi expScore, rapporto
              colonna_ricerca_D2 <- data.frame(expScore = expScore, rapporto = rapporto)
            }
            
            # costruzione righe come data.frame per il tibble finale D
            riga_ricerca_D <- data.frame(gruppo = gruppo, identificativo = identificativo, colonne_ricerca_D, colonna_ricerca_D2)
            # aggiunta righe al tibble finale D
            final_tibble_D <- bind_rows(final_tibble_D, riga_ricerca_D)
          }
        }
      }
      # scrittura degli oggetti nel file di Log
      write("prima_ricerca",out_log,append=TRUE) 
      write("seconda_ricerca",out_log,append=TRUE) 
      write("terza_ricerca",out_log,append=TRUE) 
      write("colonne_risultato_B",out_log,append=TRUE) 
      write("quarta_ricerca",out_log,append=TRUE) 
      write("quinta_ricerca",out_log,append=TRUE) 
      write("colonne_risultato_C",out_log,append=TRUE) 
      write("colonne_risultato_D",out_log,append=TRUE) 
      write("colonne_risultato_D2",out_log,append=TRUE) 
      
      # elenco degli oggetti prodotti in output
      write("\nELENCO OGGETTI OUTPUT:",out_log,append=TRUE)
      # scrittura degli oggetti
      write("final_tibble_B",out_log,append=TRUE)
      write("final_tibble_C",out_log,append=TRUE)
      write("final_tibble_D",out_log,append=TRUE)
    }
    # disconnessione da UCSC Genome Browser nell'assembly di Uomo HG38
    dbDisconnect(mydb)
    # chiusura file di Log
    close(out_log)
    
    # scrittura degli oggetti tibble B,C,D prodotti dalla funzione nel file di output file_out
    write.table(final_tibble_B, file = file_out, col.names= TRUE)
    write("\n", file = file_out, append = TRUE)
    write.table(final_tibble_C, file = file_out, append = TRUE, col.names= TRUE)
    write("\n", file = file_out, append = TRUE)
    write.table(final_tibble_D, file = file_out, append = TRUE, col.names= TRUE)
    
    # restituzione oggetti di interesse prodotti dalla funzione
    return_List <- list("final_tibble_B" = final_tibble_B, "final_tibble_C" = final_tibble_C,  "final_tibble_D" = final_tibble_D)
    return(return_List)
  } else {
    # tutti gli eventuali segnali di errore riscontrati nell'input
    error_list=character()
    write(paste0("gruppo=",toString(gruppo)," ", toString(class(gruppo))),out_log,append=TRUE)
    write(paste0("file_id=",toString(file_id)," ", toString(class(file_id))),out_log,append=TRUE)
    write(paste0("input_number=",toString(input_number)," ", toString(class(input_number))),out_log,append=TRUE)
    write(paste0("file_out=",toString(file_out)," ", toString(class(file_out))),out_log,append=TRUE)
    
    # messaggi di errore relativi agli input
    if(!file.exists(file_id)){
      error_list=c(error_list, "file_id non trovato;")
    }
    if (class(gruppo)!="character"){
      error_list=c(error_list,"gruppo deve essere di tipo character;")
    }
    if (class(file_id)!="character"){
      error_list=c(error_list,"file_id deve essere di tipo character;")
    }
    if (class(input_number)!="numeric"){
      error_list=c(error_list,"input_number deve essere di tipo numeric;")
    }
    if ((input_number)<1 || (input_number)>53){
      error_list=c(error_list,"input_number deve essere compreso fra 1 e 53;")
    }
    if (class(file_out)!="character" || length(file_out)!=1 || grepl(".+\\.txt$",file_out)==FALSE){
      error_list=c(error_list,"file_out deve essere di tipo character e avere estensione .txt;")
    }
    # scrittura nel file di Log dei messaggi di errore prodotti dalla funzione 
    write("\nERRORI:",out_log,append=TRUE)
    write(error_list,out_log,append=TRUE)
    # chiusura file di Log
    close(out_log)
    stop(paste(error_list, collapse="\n"))
  }
}
