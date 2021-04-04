########################################################################################
# banners.cmake
#
# collection of font art for Universal

####
# macro to print out a banner
macro (print_speed_banner)
    message("")
    message("___  / / /__  __ \\__  __ \\      ___  __ )__  /___    |_  ___/")
    message("__  /_/ /__  /_/ /_  /_/ /________  __  |_  / __  /| |____ \\") 
    message("_  __  / _  ____/_  _, _/_/_____/  /_/ /_  /___  ___ |___/ /") 
    message("/_/ /_/  /_/     /_/ |_|        /_____/ /_____/_/  |_/____/") 
    message("")
endmacro (print_speed_banner)

macro (print_header)
    print_speed_banner()
endmacro (print_header)

macro (print_footer)
    print_speed_banner()
endmacro (print_footer)
