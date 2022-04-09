# Written by Junyi Chen 
# 2022/2/22 Zhejiang Unv

SRC_DIR = src
VPATH = $(SRC_DIR)
CUR_DIR = $(shell pwd)
TARGET_NAME = methods
Complier = g++
Optm = -O3
OBJ = $(patsubst %.cpp,%.o,$(wildcard $(SRC_DIR)/*.cpp))
OBJ_WITHOUT_PATH = $(notdir $(OBJ))
OUT_DIR = $(CUR_DIR)/out
OBJ_DIR = $(OUT_DIR)/obj
LIB_DIR = $(OUT_DIR)/lib
OBJ_WITH_DIR = $(addprefix $(OBJ_DIR)/,$(OBJ_WITHOUT_PATH))
TARGET = $(LIB_DIR)/$(TARGET_NAME)
all: DIR_CREATE $(TARGET)

define CRT_DIR
	if [ ! -d $(1) ];then\
	mkdir -p $(1);\
	fi
endef

DIR_CREATE: 
	@$(call CRT_DIR, $(OUT_DIR))
	@$(call CRT_DIR, $(OBJ_DIR))
	@$(call CRT_DIR, $(LIB_DIR))

$(TARGET): $(OBJ_WITH_DIR)
	$(Complier) $^ -o $@
	
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(Complier) -c $(Optm) $< -o $@

.PHONY:clean
clean:
	rm -rf $(OUT_DIR)
