# Mathew McLean
# Example code for RefManager.R
# January 30, 2013

# change fname to some document on your computer and location to the directory where it is located
# location needs to end with a '/' if it is a directory e.g. 'C:/My Documents/' not 'C:/My Documents'
# location can also be a url
location <-  # 'C:\\Documents and Settings\\'
fname <-  #  'mydoc.pdf'
topic <- 'hidden Markov models'
bibEntry <- "@article{langrock2013combining,
                   title={Combining hidden Markov models for comparing the dynamics of multiple sleep electroencephalograms},
                   author={Langrock, R. and Swihart, B.J. and Caffo, B.S. and Punjabi, N.M. and Crainiceanu, C.M.},
                   journal={Statistics in Medicine},
                   year={2013},
                   publisher={Wiley Online Library}
 }"
myDB <- add.refDB(NULL,bibEntry,location,fname,topic)
print(myDB,1)
open(myDB,1)

bibEntry <- "@article{wood2012straightforward,
  title={Straightforward intermediate rank tensor product smoothing in mixed models},
  author={Wood, S.N. and Scheipl, F. and Faraway, J.J.},
  journal={Statistics and Computing},
  pages={1--20},
  year={2012},
  publisher={Springer}
}"
topic <- 'GAMs, GAMMs'
location <- 'http://link.springer.com/article/10.1007%2Fs11222-012-9314-z'
fname <- ''
myDB <- add.refDB(myDB,bibEntry,location,fname,topic)
print(myDB,2)
open(myDB,2)
myDB
search.refDB(myDB,searchterms='wood',searchfields='author')
search.refDB(myDB,c('wood','2011'),c('author','year'))
search.refDB(myDB,c('wood','2012'),c('author','year'))
search.refDB(myDB,c('caffo','2013'),c('author','year'))
getbibtex.refDB(myDB,2)
myDB <- remove.refDB(myDB,1)
attr(myDB,'size')

