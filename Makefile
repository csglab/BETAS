MAKE = make		#change this line if you are using a different GNU make software

dirBETAS = ./src/BETAS

all: MK_dir CC_BETAS RM_objectFiles

MK_dir:
	mkdir -p ./bin

CC_BETAS: $(dirBETAS)/Makefile
	$(MAKE) -C $(dirBETAS)
	
RM_objectFiles:
	rm -f $(dirBETAS)/*.o

clean:
	rm -f $(dirBETAS)/*.o
