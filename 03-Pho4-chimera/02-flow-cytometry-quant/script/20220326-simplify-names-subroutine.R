# from https://rdrr.io/cran/PTXQC/src/R/fcn_misc.R
#' Compute longest common substring of two strings.
#' 
#' Implementation is very inefficient (dynamic programming in R)
#' --> use only on small instances
#' 
#' @param s1 String one
#' @param s2 String two
#' @return String containing the longest common substring
#' 
#' @export
LCS <- function(s1, s2) 
{
  if (nchar(s1)==0 || nchar(s2)==0) return("")
  
  
  v1 <- unlist(strsplit(s1,split=""))
  v2 <- unlist(strsplit(s2,split=""))
  
  
  num <- matrix(0, nchar(s1), nchar(s2))    
  maxlen <- 0
  pstart = 0
  
  for (i in 1:nchar(s1)) {
    for (j in 1:nchar(s2)) {
      if (v1[i] == v2[j]) {
        if ((i==1) || (j==1)) { 
          num[i,j] <- 1
        } 
        else {
          num[i,j] <- 1+num[i-1,j-1]
        }
        if (num[i,j] > maxlen) {
          maxlen <- num[i,j]
          pstart = i-maxlen+1
        }
      }
    }
  }
  
  ## return the substring found
  return (substr(s1, pstart, pstart+maxlen-1))
}

#' Find longest common substring from 'n' strings.
#' 
#' Warning: greedy heuristic! This is not guaranteed to find the best solution (or any solution at all), since its done pairwise with the shortest input string as reference.
#' 
#' @param strings A vector of strings in which to search for LCS
#' @param min_LCS_length Minimum length expected. Empty string is returned if the result is shorter
#' @return longest common substring (or "" if shorter than \code{min_LCS_length})
#' 
#' @examples
#' LCSn(c("1_abcde...",
#'        "2_abcd...",
#'        "x_abc..."))  ## --> "_abc"
#' LCSn(c("16_IMU008_CISPLA_E5_R11", 
#'        "48_IMU008_CISPLA_P4_E7_R31",
#'        "60_IMU008_CISPLA_E7_R11"), 3) ## -->"_IMU008_CISPLA_"
#' LCSn(c("AAAAACBBBBB", 
#'        "AAAAADBBBBB",
#'        "AAAABBBBBEF",
#'        "AAABBBBBDGH")) ## -->  "BBBBB"
#' LCSn(c("AAAXXBBB",
#'        "BBBXXDDD",
#'        "XXAAADDD")) ## --> fails due to greedy approach; should be "XX"
#' 
#' @export
#' 
LCSn = function(strings, min_LCS_length = 0)
{
  ## abort if there is no chance of finding a suitably long substring
  if (min(nchar(strings)) < min_LCS_length) return("");
  
  if (length(strings) <= 1) return (strings)
  
  ## apply LCS to all strings, using the shortest as reference
  idx_ref = which(nchar(strings)==min(nchar(strings)))[1]
  strings_other = strings[-idx_ref]  
  candidates = unique(sapply(strings_other, LCS, strings[idx_ref]))
  candidates
  ## if only one string remains, we're done
  if (length(candidates) == 1)
  {
    if (nchar(candidates[1]) < min_LCS_length) return("") else return (candidates[1])
  }
  ## if its more, call recursively until a solution is found
  ## continue with the shortest candidates (since only they can be common to all)
  cand_short = candidates[nchar(candidates) == min(nchar(candidates))]
  solutions = unique(sapply(cand_short, function(cand) LCSn(c(cand, strings_other))))
  ## get the longest solution
  idx_sol = which(nchar(solutions)==max(nchar(solutions)))[1];
  
  return (solutions[idx_sol])
}


#'
#' Removes common substrings (infixes) in a set of strings.
#' 
#' Usually handy for plots, where condition names should be as concise as possible.
#' E.g. you do not want names like 
#' 'TK20130501_H2M1_010_IMU008_CISPLA_E3_R1.raw' and 
#' 'TK20130501_H2M1_026_IMU008_CISPLA_E7_R2.raw' 
#' but rather 'TK.._010_I.._E3_R1.raw' and
#'            'TK.._026_I.._E7_R2.raw'
#'            
#' If multiple such substrings exist, the algorithm will remove the longest first and iterate
#' a number of times (two by default) to find the second/third etc longest common substring.
#' Each substring must fulfill a minimum length requirement - if its shorter, its not considered worth removing
#' and the iteration is aborted.
#' 
#' @param strings          A vector of strings which are to be shortened
#' @param infix_iterations Number of successive rounds of substring removal
#' @param min_LCS_length   Minimum length of the longest common substring (default:7, minimum: 6)
#' @param min_out_length   Minimum length of shortest element of output (no shortening will be done which causes output to be shorter than this threshold)
#' @return A list of shortened strings, with the same length as the input                       
#'
#' @examples
#' #library(PTXQC)
#' simplifyNames(c('TK20130501_H2M1_010_IMU008_CISPLA_E3_R1.raw',
#'                 'TK20130501_H2M1_026_IMU008_CISPLA_E7_R2.raw'), infix_iterations = 2)
#' # --> "TK.._010_I.._E3_R1.raw","TK.._026_I.._E7_R2.raw"
#' 
#' try(simplifyNames(c("bla", "foo"), min_LCS_length=5))
#' # --> error, since min_LCS_length must be >=6
#' 
#' @export
#' 
simplifyNames = function(strings, infix_iterations = 2, min_LCS_length = 7, min_out_length = 7)
{
  if (min_LCS_length<6) stop( "simplifyNames(): param 'min_LCS_length' must be 6 at least.")
  
  for (it in 1:infix_iterations)
  {
    lcs = LCSn(strings, min_LCS_length=min_LCS_length)
    if (nchar(lcs)==0) return (strings) ## no infix of minimum length found
    ## replace infix with 'ab..yz' -- modified HB 2022/03/26, just remove
    strings_i = sub(paste0("(.*)", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", lcs), "(.*)"), 
                    paste0("\\1"
                           # , ifelse(substring(lcs,1,2) == "..", "", substring(lcs,1,2)) ## dont keep ".."
                           # , ".."
                           # , ifelse(substring(lcs, nchar(lcs)-1) == "..", "", substring(lcs, nchar(lcs)-1)) ## dont keep ".."
                           , "\\2"), 
                    strings)
    if (min(nchar(strings_i)) < min_out_length) return(strings) ## result got too short...
    strings = strings_i
  }
  return (strings)
}