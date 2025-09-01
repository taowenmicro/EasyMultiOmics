.onAttach <- function(libname, pkgname) {
  # 获取包版本和当前年份
  pkg_version <- utils::packageVersion("EasyMultiOmics")
  current_year <- format(Sys.Date(), "%Y")



  ascii_logo <- c(
    "  ______                __  __       _ _   _     ___              _              ",
    " |  ____|              |  \\/  |     | | | (_)   / _ \\            (_)             ",
    " | |__   __ _ ___ _   _| \\  / |_   _| | |_ _ _ | | | | _ __ ___   _   ___   ___  ",
    " |  __| / _` / __| | | | |\\/| | | | | | __| | '| | | || '_ ` _ \\ | | / __| / __| ",
    " | |___| (_| \\__ \\ |_| | |  | | |_| | | |_| | || |_| || | | | | || || (__  \\__ \\ ",
    " |______\\__,_|___/\\__, |_|  |_|\\__,_|_|\\__|_|_| \\___/ |_| |_| |_||_\\_|||___//__/ ",
    "                   __/ |                                                         ",
    "                  |___/                                                          ",
    ""
  )




  packageStartupMessage(
    paste(ascii_logo, collapse = "\n"), "\n",
    "EasyMultiOmics v", pkg_version, " (", current_year, ")\n",
    "Developed by Sheng Qirong's Microbiome Research Group\n",
    "\n",
    "Maintainers:\n",
    "  - Tao Wen <taowen@njau.edu.cn>\n",
    "  - Peng-Hao Xie <2019103106@njau.edu.cn>\n",
    "\n",
    "GitHub: https://github.com/taowenmicro/EasyMultiOmics\n",
    "Bug Reports: https://github.com/taowenmicro/EasyMultiOmics/issues\n",
    "\n",
    "Type citation(\"EasyMultiOmics\") for how to cite this package.\n",
    "Run browseVignettes(\"EasyMultiOmics\") for documentation.\n"
  )
}
