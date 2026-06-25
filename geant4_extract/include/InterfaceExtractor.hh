#pragma once

#include "DetectorAssembly.hh"

#include <map>
#include <string>

class InterfaceExtractor {

public:
    void Extract(DetectorAssembly& assembly,
                 const std::map<std::string, int>& optical_detectors,
                 double fuzzy_mm);

    void WriteInterfacesJSON(const DetectorAssembly& assembly,
                             const std::string& outDir);
};
