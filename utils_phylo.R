
# function to create nexus as string
safe_nexus = function(data, f){
  header = generate_header(data)
  body = generate_body(data)
  tail = get_tail()
  
  # safe nexus file
  cat(header, body, tail, file=f, sep = "\n")
}


# generate nexus file header
generate_header = function(data){
  header1 = "#NEXUS
BEGIN data;[eröffnet den Data-Block]
Dimensions ntax="
  header2 = " nchar="
  header3="; [Definiert die Größe des Alignments]
Format
datatype=standard
missing=?
gap=- [Definiert den Datentyp und Symbole für fehlende Daten (?) und gaps (-)]
Symbols=\"0 1 2 3 4 5 6 7 8 9\";
Matrix [hier beginnt das Alignment...]"
  header = paste(header1, nrow(data), header2, ncol(data)-1, header3, sep = "")
  header
}


normalize_name = function(x){
  n = as.character(x)
  n = paste(n, "____________", sep = "")
  n = substr(n, start = 1, stop = 10)
  n = paste(n, " ", sep = "")
  n
}


# generate data alignment string from data
generate_body = function(data){
  lines = c()
  for (r in 1:nrow(data)){
    line = normalize_name(data[r,1])
    for (c in 2:ncol(data)){
      line = paste(line, as.character(data[r,c]), sep="")
    }
    lines = append(lines, line)
  }
  body = paste(lines, sep = "\n")
}


get_tail = function(){
  tail = "; [...und hier endet es]
END; [beendet den Data-Block]"
  tail
}



# TODO make names unique




