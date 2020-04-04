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
		Tags {"Queue" = "Transparent+1" "RenderType" = "Transparent" "DisableBatching" = "True"}
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

			uniform StructuredBuffer<float4> mrl : register(t3);

			uniform sampler2D colormap;

			uniform float _Scale;

			uniform int offset;
			uniform int NumBodiesMax;
			uniform float minSize;
			uniform float maxSize;

			uniform float BHSize;
			uniform float NSSize;
			uniform float WDSize;
			uniform float NmlSize;

			struct appdata
			{
				float4 vertex : POSITION;
				UNITY_VERTEX_INPUT_INSTANCE_ID
			};

			struct v2f {
				float4 pos : SV_POSITION;
				float2 screenPos : SCREEN_POSITION;
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

  




				//scaling, if desired
				o.rad = mrl[instancePos + offset].z;
				o.teff = 0.0;
				o.special = 2;

				float rad = 1.;
				if (mrl[instancePos + offset].w > 0.) {
					rad = clamp(log10(mrl[instancePos + offset].z)*_Scale*NmlSize, minSize, maxSize);
					o.teff = TeffFromLR(mrl[instancePos + offset].w, mrl[instancePos + offset].z);
					o.special = 0;
				} else{

					if (mrl[instancePos + offset].w == -40.){ //supernovae
						o.color = float4(1., 1., 0., 1.);
						rad = 0.5;
						o.special = 1;
					}

					if (mrl[instancePos + offset].w == -30.){ //black holes
						o.color = float4(0., 1., 0., 0.5);
						rad = BHSize*_Scale;
					}

					if (mrl[instancePos + offset].w == -20.){ //neutron stars
						o.color = float4(1., 0., 1., 0.5);
						rad = NSSize*_Scale;
					}

					if (mrl[instancePos + offset].w == -10.){ //white dwarfs
						o.color = float4(0., 1., 1., 0.5);
						rad = WDSize*_Scale;
					}


				}



				// This gives us the camera's origin in 3D space (the position (0,0,0) in Camera Space)
				float4 camPos = float4(UnityObjectToViewPos(float3(0.0, 0.0, 0.0)), 1.0);  

				// extract world pivot position from object to world transform matrix
				float3 worldPos = unity_ObjectToWorld._m03_m13_m23;
 
				// apply transform scale to xy vertex positions
				float2 vertex = v.vertex.xy*rad*5.;
 
 				float2 HR = float2(log10(clamp(o.teff, 1e-5, 1e5)), log10(clamp(mrl[instancePos + offset].w, 1e-10, 1e10)));
 				HR.x = (5. - HR.x)*2. - 6.5; //scaling and repositioning -- this should be made more general and be part of the input
 				HR.y = HR.y*0.3 - 4.5; //scaling and repositioning -- this should be made more general and be part of the input

				float camDist = distance(worldPos, _WorldSpaceCameraPos);
				float4 pos = camPos + float4(HR + vertex.xy, camDist - 10.0, 0.0);
				o.pos = mul(UNITY_MATRIX_P, pos);

				float4 s = ComputeScreenPos(o.pos); 
				o.screenPos = s.xy/s.w;

				return o;
			}

			// fragment shader
			fixed4 frag(v2f i) : SV_Target
			{
				UNITY_SETUP_INSTANCE_ID(i);
				if (i.rad <= 0.) discard;

				//limit to be within the plot window.  This should also be set on input, and also the screen resolution
				if (i.screenPos.y*1080 > 300. || i.screenPos.x*1920 < 300.) discard;

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
				if (i.special >= 2 && dist >= 0.9) alpha *=1.5;

				return float4(useColor, alpha);
			}
			ENDCG
		}
	}
}