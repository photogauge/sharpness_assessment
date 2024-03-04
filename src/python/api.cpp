#include "sharpness_assessment.h"

#include <pybind11/pybind11.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

PYBIND11_MODULE(_python_api, m) {
    m.doc() = R"pbdoc(
        Python wrapper for `sharpness_assessment`.

        This information will be displayed when using `help()`:
        $ python -c "import sharpness_assessment; help(sharpness_assessment)"
    )pbdoc";

    m.def("run_sharpness_assessment", &run_sharpness_assessment, "Takes an image and assess the sharpness");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
