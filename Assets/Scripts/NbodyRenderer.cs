//https://gitlab.com/tzaeru/Nbody-simulation-tutorial-code-project
//https://www.tzaeru.com/site/compute-shader-basics-with-unity/
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class NbodyRenderer : MonoBehaviour {
	public Material material;

	public Texture colormap;
	public float minSize = 0.02f;
	public float maxSize = 0.5f;
	public float BHSize = 0.1f;
	public float NSSize = 0.1f;
	public float WDSize = 0.1f;

	Matrix4x4[][] transformList;
	int[] instances;
	int Ninstances;

	Mesh mesh;

	//This is NumBodies from the NbodyCompute C# script, I need to link these below
	private static int NumBodies;

	//sounds like it can only render 1023 at a time
	const int instance_max = 1023;

	private Transform CameraTarget;

    public void BHSizeReceiver(float val){
        BHSize = val;
    }

	// Use this for initialization
	void Start () {
		CameraTarget = GameObject.Find("CameraController").GetComponent<MouseOrbitImproved>().target;

		NbodyCompute nb = GetComponent<Transform>().gameObject.GetComponent<NbodyCompute>();
		NumBodies = nb.NumBodies;

		Ninstances = Mathf.Max(NumBodies/instance_max,1);

		transformList = new Matrix4x4[Ninstances][];
		instances = new int[Ninstances];

		MeshFilter mf = GetComponent<Transform>().gameObject.AddComponent<MeshFilter>();

		mesh = new Mesh();
		//mesh.bounds = new Bounds(Vector3.zero, Vector3.one*10000.0f);
		//mesh.bounds = new Bounds(new Vector3(23454.93f, 0.0f, 0.0f), Vector3.one*100000.0f);
		mf.mesh = mesh;

		// Create a basic quad.
		Vector3[] vertices = new Vector3[4];

		vertices[0] = new Vector3(-1, -1, 0);
		vertices[1] = new Vector3( 1, -1, 0);
		vertices[2] = new Vector3(-1,  1, 0);
		vertices[3] = new Vector3( 1,  1, 0);

		mesh.vertices = vertices;

		int[] tri = new int[6];

		tri[0] = 0;
		tri[1] = 2;
		tri[2] = 1;

		tri[3] = 2;
		tri[4] = 3;
		tri[5] = 1;

		mesh.triangles = tri;

		// Only 1023 objects can be rendered as instanced.
		// So we split the objects into sets of size 1023.
		for (int set = 0; set < Ninstances; set++)
		{
			instances[set] = instance_max;
			if (set == Ninstances - 1)
			{
				instances[set] = Mathf.Min(NumBodies, NumBodies % instance_max); //to account for <1023 objects
			}

			transformList[set] = new Matrix4x4[instances[set]];

			for (int i = 0; i < instances[set]; i++)
			{
				Matrix4x4 matrix = new Matrix4x4();
				//matrix.SetTRS(Vector3.zero, Quaternion.Euler(Vector3.zero), Vector3.one);
				//THIS SETTING IN VERY IMPORTANT.  Instancing will only be active (and therefore the object will only be drawn) when the first Vector3 is within the camera's view.  For instance, if it is set to Vector3.zero, the Nbody particles will only be drawn when the scene's origin is in view of the camera!
				matrix.SetTRS(CameraTarget.position, Quaternion.Euler(Vector3.zero), Vector3.one);
				transformList[set][i] = matrix;
			}
		}
	}
	
	// Update is called once per frame
	void Update () {

		for (int set = 0; set < Ninstances; set++){
			MaterialPropertyBlock mpb = new MaterialPropertyBlock();
			mpb.SetInt("offset", set*instance_max);
			mpb.SetFloat("minSize", minSize);
			mpb.SetFloat("maxSize", maxSize);
			mpb.SetFloat("BHSize", BHSize);
			mpb.SetFloat("NSSize", NSSize);
			mpb.SetFloat("WDSize", WDSize);
			mpb.SetTexture("colormap", colormap);

			for (int i = 0; i < instances[set]; i++)
			{
				Matrix4x4 matrix = new Matrix4x4();
				//matrix.SetTRS(Vector3.zero, Quaternion.Euler(Vector3.zero), Vector3.one);
				//THIS SETTING IN VERY IMPORTANT.  Instancing will only be active (and therefore the object will only be drawn) when the first Vector3 is within the camera's view.  For instance, if it is set to Vector3.zero, the Nbody particles will only be drawn when the scene's origin is in view of the camera!
				matrix.SetTRS(CameraTarget.position, Quaternion.Euler(Vector3.zero), Vector3.one);
				transformList[set][i] = matrix;
			}


			Graphics.DrawMeshInstanced(mesh, 0, material, transformList[set], instances[set], mpb);
		}
	}
}
