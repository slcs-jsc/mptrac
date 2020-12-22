#define RED 0xFFFF0000
#define BLUE 0xFF0000FF
#define GREEN 0xFF008000
#define YELLOW 0xFFFFFF00
#define CYAN 0xFF00FFFF
#define MAGENTA 0xFFFF00FF
#define GRAY 0xFF808080
#define PURPLE 0xFF800080

// Macro for calling nvtxRangePushEx
#define RANGE_PUSH(range_title,range_color) { \
    nvtxEventAttributes_t eventAttrib = {0}; \
    eventAttrib.version = NVTX_VERSION; \
    eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE; \
    eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;\
    eventAttrib.colorType = NVTX_COLOR_ARGB; \
    eventAttrib.color = range_color; \
    eventAttrib.message.ascii = range_title; \
nvtxRangePushEx(&eventAttrib); \
}

// Macro for calling nvtxRangePop
#define RANGE_POP {\
    nvtxRangePop();\
}
