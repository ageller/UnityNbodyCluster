using UnityEngine;
using UnityEngine.UI;
using System.Collections;
using UnityEngine.SceneManagement;

public class ButtonController : MonoBehaviour {

	private bool started = false;
	private bool paused = false;

	void Start(){}

	public void StartButtonClicked(){
		if (started){
			if (paused){
				//Debug.Log("playing...");
				paused = false;
				Time.timeScale = 1;
				GameObject.Find("StartButton").GetComponentInChildren<Text>().text = "Pause";

				//do something here

			} else {
				//Debug.Log("pausing...");
				paused = true;
				Time.timeScale = 0;
				GameObject.Find("StartButton").GetComponentInChildren<Text>().text = "Play";

				//do something here
			}
		} else {
			//Debug.Log("starting game...");
			started = true;
			Time.timeScale = 1;
			GameObject.Find("StartButton").GetComponentInChildren<Text>().text = "Pause";

			//do something here

		}
	}

	public void RestartButtonClicked(){
		SceneManager.LoadScene("NbodySimulator");
	}
}
