using UnityEngine;
using UnityEngine.UI;
using System.Collections;

public class BHSizeSliderController : MonoBehaviour {
	
	private NbodyRenderer renderer;
	private Slider slider;

	private float value = 1.0f;
	

	void Awake(){
		renderer = GameObject.Find("NbodyCompute").GetComponent<NbodyRenderer>();
		slider = GetComponent<Slider>();
		value = renderer.BHSize;
		slider.value = renderer.BHSize;

	}

	void Update(){

		if (slider.value != value) {
			value = slider.value;
			renderer.SendMessage("BHSizeReceiver", value);

		}

	}

}
