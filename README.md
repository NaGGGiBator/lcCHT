# lcCHT

New solver for Openfoam

## **Include**

functions
{
    libs
    (
	"fieldFunctionObjects"
	"lc_cht_lib.so"
    );

    lcCHT
    {
        type            	lcCHT;
	SETTING
    }
}


## **Settings**
**Important parametrs:**
*region* — The region where to look for the patch
*patchName* — Coupling interface
*time_step_fluid* — Time step in the fluid zone
*time_step_solid* — Time step in the solid zone

**Optional:**

### **Update 18.11.2025**

