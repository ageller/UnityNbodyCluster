using UnityEngine;
using System.Collections;

public class GUITextDisplay : MonoBehaviour
{

	private MouseOrbitImproved mainCamera;
	float deltaTime = 0.0f;
	float time = 0.0f;

	void Start(){
		//camera = GameObject.Find("CameraController").GetComponent<MouseOrbitImproved>();
		mainCamera = GameObject.Find("MainCamera").GetComponent<MouseOrbitImproved>();

	}

	void Update(){
		deltaTime += (Time.unscaledDeltaTime - deltaTime) * 0.1f;
		time += Time.deltaTime;

	}
	void OnGUI(){
		int w = Screen.width, h = Screen.height;
		int fs = h*2/100;

		GUIStyle style = new GUIStyle();
		Rect rect = new Rect(0f,0f,0f,0f);
		string text = "";

		float y = 10f;
		float x = 10f;

		//Time
		rect = new Rect(x, y, w, fs);
		style.alignment = TextAnchor.UpperLeft;
		style.fontSize = fs;
		style.normal.textColor = new Color(1.0f, 1.0f, 1.0f, 1.0f);
		text = string.Format("{0:0.0} Myr", time);
		GUI.Label(rect, text, style);

		//Distance
		y += 1.2f*fs;
		rect = new Rect(x, y, w, fs);
		style.alignment = TextAnchor.UpperLeft;
		style.fontSize = fs;
		style.normal.textColor = new Color(1.0f, 1.0f, 1.0f, 1.0f);
		text = string.Format("{0:0.0} pc", mainCamera.distance);
		GUI.Label(rect, text, style);

		//FPS
		y+= 1.2f*fs;
		rect = new Rect(x, y, w, fs);
		style.alignment = TextAnchor.UpperLeft;
		style.fontSize = fs;
		style.normal.textColor = new Color(1.0f, 1.0f, 1.0f, 1.0f);
		float msec = deltaTime*1000.0f;
		float fps = 1.0f/deltaTime;
		text = string.Format("{0:0.0} ms ({1:0.} fps)", msec, fps);
		GUI.Label(rect, text, style);

	}
}