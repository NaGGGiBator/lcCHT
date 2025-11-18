# lcCHT

New solver for Openfoam

## **Include**

```C++
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
```

## **Settings**
**Important parametrs:** 				<br>
*region* — The region where to look for the patch	<br>
*patchName* — Coupling interface			<br>
*time_step_fluid* — Time step in the fluid zone		<br>
*time_step_solid* — Time step in the solid zone		<br>

**Optional:**

### **Update 18.11.2025**

