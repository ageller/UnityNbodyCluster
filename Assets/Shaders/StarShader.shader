//https://gitlab.com/tzaeru/nbody-simulation-tutorial-code-project
//https://www.tzaeru.com/site/compute-shader-basics-with-unity/
Shader "Custom/Star"
{

    Properties
    {
        _Scale ("Scale", Float) = 1.0
    }
	SubShader
	{
		//added disable batching to stop flickering
		Tags {"Queue" = "Transparent" "RenderType" = "Transparent" "DisableBatching" = "True"}
		LOD 100

		Cull Off
		ZWrite Off
		AlphaTest Off
		Blend SrcAlpha OneMinusSrcAlpha

		Pass
		{
			CGPROGRAM
			
			#pragma target 5.0

			#pragma vertex vert
			#pragma fragment frag

			#pragma multi_compile_instancing

			#include "UnityCG.cginc"

			//I think this is where the position is pulled from the compute shader
			//uniform StructuredBuffer<float3> position : register(t1);
			uniform StructuredBuffer<float3> positionCoM : register(t8);
			uniform StructuredBuffer<float> radius : register(t4);
			uniform StructuredBuffer<float> luminosity : register(t5);

			uniform sampler2D colormap;

			uniform int offset;
			uniform float _Scale;
			uniform float minSize;
			uniform float maxSize;

			uniform float BHSize;
			uniform float NSSize;
			uniform float WDSize;
			uniform float NmlSize;

			//in units of parsecs
			//#define RSun 2.2546101516841093e-08

			struct appdata
			{
				float4 vertex : POSITION;
				UNITY_VERTEX_INPUT_INSTANCE_ID
			};

			struct v2f {
				float4 pos : SV_POSITION;
				float4 vertex : VERTEX;
				float rad: RADIUS;
				float teff: TEFF;
				float4 color: COLOR;
				int special: SPECIAL;
				UNITY_VERTEX_INPUT_INSTANCE_ID
			};



			float TeffFromLR(float L, float R){
			//use stellar radius and stellar luminosity to get the star's effective temperature
				if (L > 0.){
					float logTeff = 3.762f + 0.25f*log10(L) - 0.5f*log10(R);
					return pow(10.0f, logTeff);
				} else return L;
			}

			// Vertex shader
			v2f vert(appdata v)
			{
				v2f o;
				UNITY_SETUP_INSTANCE_ID(v);
				UNITY_TRANSFER_INSTANCE_ID(v, o);
				
				//add the actual position from the instance, keeping in view space
				int instancePos = 0;
#ifdef UNITY_INSTANCING_ENABLED
				instancePos = unity_InstanceID;
#endif

				o.vertex = v.vertex;

				//https://forum.unity.com/threads/programmer-friendly-billboard-shader.724946/
				// This gives us the camera's origin in 3D space (the position (0,0,0) in Camera Space)
				// Note: UnityObjectToViewPos(pos) is equivalent to mul(UNITY_MATRIX_MV, float4(pos, 1.0)).xyz,
				float4 camPos = float4(UnityObjectToViewPos(float3(0.0, 0.0, 0.0)), 1.0);    

				// Since w is 0.0, in homogeneous coordinates this represents a vector direction instead of a position
				float4 viewDir = float4(v.vertex.x, v.vertex.y, 0.0, 0.0);  

				//scaling, if desired
				float rad = 1.;
				if (luminosity[instancePos + offset] > 0.) {
					rad = clamp(log(radius[instancePos + offset])*_Scale, minSize, maxSize)*NmlSize;
					o.teff = TeffFromLR(luminosity[instancePos + offset], radius[instancePos + offset]);
					o.special = 0;
				} else{

					o.special = 2;
					o.teff = 0.;

					if (luminosity[instancePos + offset] == -40.){ //supernovae
						o.color = float4(1., 1., 0., 1.);
						rad = 3.;
						o.special = 1;
					}

					if (luminosity[instancePos + offset] == -30.){ //black holes
						o.color = float4(0., 1., 0., 0.5);
						rad = BHSize*_Scale;
					}

					if (luminosity[instancePos + offset] == -20.){ //neutron stars
						o.color = float4(1., 0., 1., 0.5);
						rad = NSSize*_Scale;
					}

					if (luminosity[instancePos + offset] == -10.){ //white dwarfs
						o.color = float4(0., 1., 1., 0.5);
						rad = WDSize*_Scale;
					}


				}


				o.rad = rad;
				
				//float rad = 1.0;
				float4 scale = float4(rad, rad, 1.0, 1.0); 
				
				// Add the camera position and direction, 
				float4 pos = camPos + viewDir*scale; 

				pos += float4(UnityObjectToViewPos(positionCoM[instancePos + offset]), 0.0f);
				//pos += float4(UnityObjectToViewPos(position[instancePos + offset]), 0.0f);

				//then multiply by UNITY_MATRIX_P to get the new projected vertex position
				o.pos = mul(UNITY_MATRIX_P, pos);

				return o;
			}

			// fragment shader
			fixed4 frag(v2f i) : SV_Target
			{
				UNITY_SETUP_INSTANCE_ID(i);
				if (i.rad <= 0.) discard;

				float dist = distance(i.vertex, float4(0.0, 0.0, 0.0, 0.0))/sqrt(2.0);

				if (dist >= 1.) discard;

				float3 useColor = i.color.rgb;
				float alpha = i.color.a;

				if (i.special < 1){
					float scaledTeff = clamp( ((i.teff - 1000.)/19000.), 0.01, 0.99);
					float4 color = tex2D(colormap, float2(scaledTeff, 0.5));
					float brightness = 1. - dist;
					useColor = float3(lerp(float3(1.0, 1.0, 1.0), color.rgb, clamp(1.1 - brightness, 0.0, 1.0)));
				}
				if (i.special < 2){
					alpha = 1. - pow(dist, 3.);
				}

				return float4(useColor, alpha);
			}
			ENDCG
		}
	}
}