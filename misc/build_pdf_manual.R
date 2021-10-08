pkg_path = find.package('MiaoCom')
system(paste(shQuote(file.path(R.home("bin"), "R")),
             "CMD",
             "Rd2pdf",
             shQuote('D:/Program Files/R/R-4.1.0/library/MiaoCom')))
