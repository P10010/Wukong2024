#include "../include/App.h"
#include "../include/PBD.h"

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <path to obj file>" << std::endl;
        return 1;
    }

    PBD pbd;

    pbd.initializeFromFile(argv[1]);

    App<PBD> app(pbd);
    app.initializeScene();
    app.run();

    return 0;
}
