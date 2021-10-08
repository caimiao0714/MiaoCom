pkg_path = find.package('MiaoCom')
system(paste(shQuote(file.path(R.home("bin"), "R")),
             "CMD",
             "Rd2pdf",
             shQuote(pkg_path)))
