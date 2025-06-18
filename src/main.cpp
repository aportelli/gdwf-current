#include <Grid/Grid.h>

using namespace Grid;

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  std::cout << GridLogMessage << "Grid initialized successfully." << std::endl;

  // Your application logic here

  Grid_finalize();

  std::cout << GridLogMessage << "Grid finalized successfully." << std::endl;

  return EXIT_SUCCESS;
}
