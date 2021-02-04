import React, { Component } from 'react'

export class Form extends Component
{
    state = {
        name: '',
        example_int: '',
    }

    onChange = e => this.setState({ [e.target.name]: e.target.value });
    
    onSubmit = e => {
        e.preventDefault();
        console.log(this.state)
    }

    render()
    {
        const { name, example_int} = this.state;
        return (
            <div className="card card-body mt-4 mb-4">
                <h2>Example Form</h2>
                <form onSubmit={this.onSubmit}>
                    <div className="form-group">
                        <label>Fedid</label>
                            <input
                                className="form-control"
                                type="text"
                                name="name"
                                onChange={this.onChange}
                                value={name}
                            />
                    </div>

                    <div className="form-group">
                        <label>Example Integer</label>
                        <select className="form-control" id="di_sel" name="decision_int" onChange={this.onChange} value={example_int}>
                            <option>1</option>
                            <option>2</option>
                            <option>3</option>
                            <option>4</option>
                        </select>
                    </div>
                    <div className="form-group">
                        <button type="submit" className="btn btn-primary">Submit</button>
                    </div>
                </form>
            </div>
        )
    }
}

export default Form