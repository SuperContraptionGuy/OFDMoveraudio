BUILD_DIR := build

CC := gcc

SRCS := ofdm.c ofdmDecoder.c
OBJS := $(addprefix $(BUILD_DIR)/,$(SRCS:.c=.o))
DEPS = $(OBJS:%.o=%.d)

DEPFLAGS = -MT "$@" -MMD -MP -MF "$(BUILD_DIR)/$*.d"

CFLAGS := \
	-Wall \
	-Wextra \
	-Wstrict-prototypes \
	-Wconversion \
	-Wsign-conversion \
	-Wformat-security \
	-Wstack-protector \
	-Wno-unused-function \
	-flto \
	-g3 \
	-std=gnu18 \
	-march=native \
	-mtune=native \
	-fstack-protector-all \
	--param ssp-buffer-size=4 \
	-fstack-clash-protection \
	-fPIE \
	-fPIC \
	-ftrapv \
	-D_FORTIFY_SOURCE=2 \
	-fcf-protection=full \
	-fsanitize=address \
	# -O3

LDFLAGS := \
	-Wl,-z,relro,-z,now \
	-lm \
	-lfftw3

.PHONY: all ofdm ofdmDecoder clean
all: ofdm ofdmDecoder

ofdm: $(BUILD_DIR)/ofdm

ofdmDecoder: $(BUILD_DIR)/ofdmDecoder

$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	@echo [ CC ] $@
	@$(CC) -x c $(CFLAGS) $(DEPFLAGS) -c $< -o $@

$(BUILD_DIR)/ofdm: $(BUILD_DIR)/ofdm.o $(BUILD_DIR)/utilities.o
	@echo [ LD ] $@
	@$(CC) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)/ofdmDecoder: $(BUILD_DIR)/ofdmDecoder.o $(BUILD_DIR)/utilities.o
	@echo [ LD ] $@
	@$(CC) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(BUILD_DIR)

$(DEPS): ;

include $(DEPS)
