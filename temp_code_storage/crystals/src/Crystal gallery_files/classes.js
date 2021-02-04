/* an element that can be shown or hidden, together with controls used to show or hide it */
class Hideable {
	constructor(container, target, targetDisplay, hiderElementSelector, showerElementSelector, controlsDisplay, targetHide = null ){
		this.container = container;								//common parent element of the target and the controls (selector)									
		this.target = target;									//element to be shown or hidden (selector)
		this.targetDisplay = targetDisplay;						//CSS display attribute or class of target when shown (string)
		this.hiderElementSelector = hiderElementSelector; 		//control used to hide target (selector)
		this.showerElementSelector = showerElementSelector;		//control used to show target (selector)
		this.controlsDisplay = controlsDisplay;					//CSS display attribute of controls when shown (string)
		this.targetHide = targetHide;							////class of target when hidden (string)
	}
	
	//use the provided selector to access the target element
	findTarget() {
		this.targetElement = this.container.querySelector(this.target);
	}
	
	//use provided selectors to access the controls
	findControls() {
		this.hiderElement = this.container.querySelector(this.hiderElementSelector);
		this.showerElement = this.container.querySelector(this.showerElementSelector);
	}

	showOneTargetElement(element) {
		//show by changing class
		if (this.targetDisplay.includes('.')) {
			element.className = this.targetDisplay.slice(1, this.targetHide.length);
		}
		//show by changing display
		else {
			element.style.display = this.targetDisplay;
		}
	}
	
	hideOneTargetElement(element) {
		//hide by changing class
		if (this.targetHide !== null) {
			element.className = this.targetHide.slice(1, this.targetHide.length);
		}
		//hide by changing display
		else {
			element.style.display = 'none';
		}
	}
	
	showTarget(){
		this.showOneTargetElement(this.targetElement);
	}
	
	hideTarget() {
		this.hideOneTargetElement(this.targetElement);
	}
	//show or hide target by clicking on controls; swap display attr of controls
	addListeners() {
		this.findTarget();
		this.findControls();
		
		this.showerElement.addEventListener('click', () => {
			this.showTarget();
			this.hiderElement.style.display = this.controlsDisplay
			this.showerElement.style.display = 'none';
			
		})
		
		this.hiderElement.addEventListener('click', () => {
			this.hideTarget();
			this.hiderElement.style.display = 'none'; 
			this.showerElement.style.display = this.controlsDisplay;
		})		
	}
}

//ShowAndHideGroup, but the target is an array of elements
class HideableArray extends Hideable {
	
	//for elements that have their individual show/hide controls
	set individualControls(array){
		this.icShowerSelector = array[0];
		this.icHiderSelector = array[1];
		this.icDisplayType = array[2];
	}
	
	findTarget() {
		this.target = this.container.querySelectorAll(this.target);
		};
	
	showTarget() {		
		this.target.forEach(element => {
			this.showOneTargetElement(element);
			//swap individual controls if present
			if (this.icShowerSelector !== undefined){
				element.parentElement.querySelector(this.icShowerSelector).style.display = 'none';
				element.parentElement.querySelector(this.icHiderSelector).style.display = this.icDisplayType;
			}
		})
	}
	
	hideTarget() {
		this.target.forEach(element => {
			this.hideOneTargetElement(element);
			//swap individual controls if present
			if (this.icShowerSelector !== undefined){
				element.parentElement.querySelector(this.icShowerSelector).style.display = this.icDisplayType; 
				element.parentElement.querySelector(this.icHiderSelector).style.display = 'none';
			}
		})
	}
}

class MassInputChecker{
	
	constructor(container, checkablesSelector, checkerSelector, uncheckerSelector=null){
		this.container = container;	
		this.checkablesSelector = checkablesSelector;
		this.checkerSelector = checkerSelector;
		this.uncheckerSelector = uncheckerSelector;
	}
	
	findElements(){
		this.checkablesArray = this.container.querySelectorAll(this.checkablesSelector);
		this.checker = this.container.querySelector(this.checkerSelector);
		if (this.uncheckerSelector !== null){
			this.unchecker = this.container.querySelector(this.uncheckerSelector);
		}
		
	}
	
	checkAll(){
		this.checker.addEventListener('click', () => {
			this.checkablesArray.forEach(input => {
				if (input.checked !== true) {
				input.click();
				}
			})
		})
	}
	
	uncheckAll(){
		this.unchecker.addEventListener('click', () => {
			this.checkablesArray.forEach(input => {
				if (input.checked === true) {
				input.click();
				}
			})
		})
	}
	
	addListeners() {
		this.findElements();
		this.checkAll();
		if (this.uncheckerSelector !== null){
			uncheckAll();
		}
	}	
}
