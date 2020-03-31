using UnityEngine;
using UnityEngine.UI;
using System.Collections;
using UnityEngine.SceneManagement;

public class ButtonController : MonoBehaviour {

	public bool started = false;
	public bool paused = true;

	void Start(){
		Time.timeScale = 0;
	}

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
			paused = false;
			Time.timeScale = 1;
			GameObject.Find("StartButton").GetComponentInChildren<Text>().text = "Pause";

			//do something here
			//disable the initial condition sliders
			GameObject.Find("vScaleSlider").GetComponentInChildren<Slider>().interactable = false;
			GameObject.Find("RhmSlider").GetComponentInChildren<Slider>().interactable = false;
			GameObject.Find("NSlider").GetComponentInChildren<Slider>().interactable = false;
		}
	}

	public void RestartButtonClicked(){
		SceneManager.LoadScene("NbodyStarCluster");
	}
}
