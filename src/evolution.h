
#include <inttypes.h>
#include <iostream>
#include <list>
#include <vector>

#pragma pack(push, 1)



// Class to register events and write them to file
class mutation{
private:
    FILE * logfile;
    char * write_buffer;
    uint64_t write_index;
    uint64_t event_count;
    char tmp[256];

public:
    ee_log(FILE * logfile);
    void register_event(Event e, void * event_data);    
    ~ee_log();

private:
    void write(const char * data);

};

