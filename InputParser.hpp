#pragma once
#include <string>
#include <map>

namespace wf_solver {
class InputParser
{
public:
    typedef std::map<std::string, std::string> container_type;
    // Constructor
    InputParser(std::string &input_filename): input_filename_(input_filename) 
    {   this->parse();    }
    bool parse();
    void clear();
    const container_type &data(){   return this->keyword_map_; };
private:
    std::string input_filename_;
    container_type keyword_map_;
};

};
