/* 
 * File:   main.cpp
 * Author: sigve
 *
 * Created on December 26, 2012, 3:18 PM
 */

#include <cstdlib>
#include "src/MainApplication.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    MainApplication *app;

    if (argc == 2)
        app = new MainApplication(&argc, &argv, argv[1]);
    else
        app = new MainApplication(&argc, &argv, "../config.cfg");

    app->runConfiguration();
    app->finalize();

    // A little output to help notice when the program is finished executing
    //    int ret = system("zenity --info --text='Completed: Qunatum Dynamics computaion using the Standard Method(SM)'");
    //    (void)ret;
    return 0;
}

