from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps, cmake_layout


class StructureIdentificationConan(ConanFile):
    name = "structure-identification"
    version = "1.0.0"
    package_type = "static-library"
    license = "MIT"
    settings = "os", "arch", "compiler", "build_type"
    requires = (
        "coretoolkit/1.0.0",
        "spdlog/1.14.1",
        "fmt/10.2.1",
        "nlohmann_json/3.11.3",
    )
    exports_sources = "CMakeLists.txt", "include/*", "src/*"

    def layout(self):
        cmake_layout(self)

    def generate(self):
        CMakeToolchain(self).generate()
        CMakeDeps(self).generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.set_property(
            "cmake_target_name", "structure-identification::structure-identification"
        )
        self.cpp_info.libs = ["structure-identification"]
