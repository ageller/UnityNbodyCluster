//https://github.com/imifos/unity3d-workspace/tree/master/Volumetric_Galaxy
//I edited and cleaned this up so that the galaxy is centered at the origin, and scaled properly

using System.Collections.Generic;
using UnityEngine;

/*
Generates a volumentic galaxy based on the passed bitmap, using the partcile system.
Developed by Natasha CARL, @Imifos
License: Public domain, use at your own risk.
*/

public class GalaxyCloudForge : MonoBehaviour
{

	public float scale = 10.0f;
	public int stepsize = 10;

	// Assigned galaxy bitmap in Inspector
	// The galaxy bitmap used in this demo is from my favorite Space Sandbox
	// Elite Dangerous, (c) Frontier Development
	public Texture2D galaxyBitmap;

	//galaxy image material
	public Material galaxyMat;

	// Particles added to the particle system of the current transform.
	// It should use a custom material set to the default (Unity3D) Particle bitmap
	// but the Built-in "Additive (Soft)" shader.
	private List<ParticleSystem.Particle> particles = new List<ParticleSystem.Particle>();

	private float textureWidth;
	private float textureHeight;
	private float textureCenterX;
	private float textureCenterY;

	/*
	 * Generates a galaxy as cloud of particles, based on a template bitmap.
	 */
	private void BuildGalaxy()
	{
		textureWidth = galaxyBitmap.width;
		textureHeight = galaxyBitmap.height;
		textureCenterX = textureWidth/2;
		textureCenterY = textureHeight/2;

		particles.Clear(); // In case we call the method in a loop for debugging
        //int attractor = 0;

		// Generate large clouds
		// ---------------------
		for (float y = 0; y < textureHeight; y += stepsize)
			for (float x = 0; x < textureWidth; x += stepsize)
			{
				Color c = galaxyBitmap.GetPixel((int)x, (int)y);

				if (c.g < 0.01)
					continue;

				float px = x - textureWidth/2.0f;
				float pz = y - textureHeight/2.0f;

				px += Random.Range(-0.001f*textureWidth, 0.001f*textureWidth);
				pz += Random.Range(-0.002f*textureWidth, 0.002f*textureWidth);

				// 0 at the borders, 1 on one of the center axes. It forms a 'square',
				// not a circle (around the center), so it's really a huge approximation.
				// Used (below) to generate larger clouds in the center, but brighter clouds 
				// outside to galactic arms better visible.
				float axisDistanceFactor = 1 - (Mathf.Abs(x - textureCenterX)/textureWidth +
												Mathf.Abs(y - textureCenterY)/textureHeight);

				// Small emphasis on "blue only" sectors in terms of luminosity
				float blueFactor = 0f;
				if (c.r < 0.2 && c.g < 0.3 && c.g > 0.01)
					blueFactor = (1 - c.b) / 4;

				ParticleSystem.Particle p = new ParticleSystem.Particle
				{
					position = new Vector3(px, 0f, pz),
					startSize = Random.Range(100f*axisDistanceFactor, 140f*axisDistanceFactor),
					startColor = new Color(c.r, c.g, c.b, 0.05f + (1.0f - axisDistanceFactor) * 0.3f + blueFactor)
				};
				p.position *= scale;
				p.startSize *= scale*stepsize/10.0f;

				particles.Add(p);


			// 	// Generate bright small clouds and 'fake' stars
			// 	// ---------------------------------------------
			// 	px = x - textureWidth/2.0f;
			// 	pz = y - textureHeight/2.0f;
			// 	float py = Random.Range(0.8f*axisDistanceFactor, 1.2f*axisDistanceFactor);

			// 	// Clouds
			// 	// Randomised on bitmap points which have a given color level
			// 	// Note: c.g = grayscale value = brightness, which is not useful here
			// 	// as this would be centered in the bright galaxy center.
			// 	float rnd = Random.Range(0, Mathf.Clamp(20 - attractor, 1, 20));
			// 	if ((c.r > 0.2f || c.g > 0.2f || c.b > 0.2f) && rnd < 2f)
			// 	{
			// 		// Increase chance to build groups of clouds, up to 3 times
			// 		attractor += 7;
			// 		if (attractor > 25)
			// 			attractor = 0;

			// 		blueFactor = 0f;
			// 		if (c.r < 0.2 && c.g < 0.5)
			// 			blueFactor = (1 - c.b) / 3;

			// 		ParticleSystem.Particle p2 = new ParticleSystem.Particle
			// 		{
			// 			position = new Vector3(px, py, pz),
			// 			startSize = Random.Range(40f*axisDistanceFactor, 80f*axisDistanceFactor),
			// 			startColor = new Color(c.r, c.g, c.b, Random.Range(0.1f, 0.3f * axisDistanceFactor + blueFactor))
			// 		};
			// 		p2.position *= scale;
			// 		p2.startSize *= scale;

			// 		particles.Add(p2);
			// 	}
			// 	else
			// 	{
			// 		// No new cloud generated, reset random() probability
			// 		attractor = 0;
			// 	}

			// 	// Stars
			// 	rnd = Random.Range(0, 50);
			// 	if ((c.r > 0.4f || c.g > 0.2f || c.b > 0.2f) && rnd < 2.5)
			// 	{
			// 		py = Random.Range(-0.2f*axisDistanceFactor, 0.2f*axisDistanceFactor);
			// 		ParticleSystem.Particle p3 = new ParticleSystem.Particle
			// 		{
			// 			position = new Vector3(px, py, pz),
			// 			startSize = Random.Range(10f, 20f),
			// 			startColor = new Color(1, 1, 1, Random.Range(0.4f, 0.7f))
			// 		};
			// 		p3.position *= scale;
			// 		p3.startSize *= scale;
			// 		//particles.Add(p3);
			// 	}

			}


		// And finally: all particles into the scene
		GetComponent<ParticleSystem>().SetParticles(particles.ToArray(), particles.Count);

		//also add the image on a plance
		GameObject my_plane = GameObject.CreatePrimitive(PrimitiveType.Plane);
		my_plane.transform.localScale = new Vector3(3000f,0f,3000f); //not sure what to use here
		my_plane.GetComponent<Renderer>().material = galaxyMat;
		my_plane.GetComponent<Renderer>().material.mainTexture = galaxyBitmap;



		Debug.Log("Galaxy - Added Particles:" + particles.Count);
	}

	// Unity3D call back at scene start
	void Start()
	{
		BuildGalaxy();
	}

}