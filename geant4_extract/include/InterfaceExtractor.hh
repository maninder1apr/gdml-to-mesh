#pragma once

#include "DetectorAssembly.hh"
#include <string>

class InterfaceExtractor {

public:
    void Extract(DetectorAssembly& assembly);
    void WriteInterfacesJSON(const DetectorAssembly& assembly, const std::string& outDir);
};
