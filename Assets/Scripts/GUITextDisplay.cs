using UnityEngine;
using System.Collections;

public class GUITextDisplay : MonoBehaviour
{

	private MouseOrbitImproved camera;
	float deltaTime = 0.0f;
	float time = 0.0f;

	void Start(){
		camera = GameObject.Find("CameraController").GetComponent<MouseOrbitImproved>();
	}

	void Update(){
		deltaTime += (Time.unscaledDeltaTime - deltaTime) * 0.1f;
		time += Time.deltaTime;

	}
	void OnGUI(){
		int w = Screen.width, h = Screen.height;
		int fs = h*2/100;

		GUIStyle style = new GUIStyle();

		//FPS
		Rect rect = new Rect(0, 0, w, h * 2 / 100);
		style.alignment = TextAnchor.UpperLeft;
		style.fontSize = h * 2 / 100;
		style.normal.textColor = new Color(1.0f, 1.0f, 1.0f, 1.0f);
		float msec = deltaTime * 1000.0f;
		float fps = 1.0f / deltaTime;
		string text = string.Format("{0:0.0} ms ({1:0.} fps)", msec, fps);
		GUI.Label(rect, text, style);

		//Time
		rect = new Rect(0, 1.1f*fs, w, fs);
		style.alignment = TextAnchor.UpperLeft;
		style.fontSize = fs;
		style.normal.textColor = new Color(1.0f, 1.0f, 1.0f, 1.0f);
		text = string.Format("{0:0.0} Myr", time);
		GUI.Label(rect, text, style);

		//Distance
		rect = new Rect(0, 2.2f*fs, w, fs);
		style.alignment = TextAnchor.UpperLeft;
		style.fontSize = fs;
		style.normal.textColor = new Color(1.0f, 1.0f, 1.0f, 1.0f);
		text = string.Format("{0:0.0} pc", camera.distance);
		GUI.Label(rect, text, style);


	}
}