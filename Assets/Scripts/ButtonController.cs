using UnityEngine;
using UnityEngine.UI;
using System.Collections;
using UnityEngine.SceneManagement;

public class ButtonController : MonoBehaviour {

	public bool started = false;
	public bool paused = true;

	void Start(){
		Time.timeScale = 0;

		GameObject.Find("CameraCenter").GetComponent<TrailController>().enabled = false;
		TrailRenderer tr = GameObject.Find("CameraCenter").GetComponent<TrailRenderer>();
		tr.enabled = false;
		tr.widthMultiplier = 0.0f;
		tr.time = 0.0f;
		tr.emitting = false;
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
			GameObject.Find("FractalSlider").GetComponentInChildren<Slider>().interactable = false;
			GameObject.Find("RgcSlider").GetComponentInChildren<Slider>().interactable = false;

			TrailRenderer tr = GameObject.Find("CameraCenter").GetComponent<TrailRenderer>();
			GameObject.Find("CameraCenter").GetComponent<TrailController>().enabled = true;


		}
	}

	public void RestartButtonClicked(){
		SceneManager.LoadScene("NbodyStarCluster");
	}
}
