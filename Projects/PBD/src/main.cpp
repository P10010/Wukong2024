#include "../include/App.h"
#include "../include/PBD.h"

int main()
{

    PBD pbd;
    // todo this is temporary, it's probably better to copy the obj in the PBD
    // folder if we keep using it
    pbd.initializeFromFile("../../../Projects/DiscreteShell/data/grid.obj");

    App<PBD> app(pbd);
    app.initializeScene();
    app.run();

    return 0;
}
