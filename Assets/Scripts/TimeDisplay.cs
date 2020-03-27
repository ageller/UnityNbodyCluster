using UnityEngine;
using System.Collections;

public class TimeDisplay : MonoBehaviour
{
	float time = 0.0f;
	void Update()
	{
		time += Time.deltaTime;
	}

	void OnGUI()
	{
		int w = Screen.width, h = Screen.height;

		GUIStyle style = new GUIStyle();

		int fs = h*2/100;
		Rect rect = new Rect(0, 1.1f*fs, w, fs);
		style.alignment = TextAnchor.UpperLeft;
		style.fontSize = fs;
		style.normal.textColor = new Color(1.0f, 1.0f, 1.0f, 1.0f);
		string text = string.Format("{0:0.0} Myr", time);
		GUI.Label(rect, text, style);
	}
}