import React, { Component, Fragment } from 'react';
import { connect } from 'react-redux';
import PropTypes from 'prop-types';
import { getExample } from '../../actions/example';

/*
This is assuming that the model table is just a table with one column called smiles
Will obviously need updating!!
*/

export class ExampleTable extends Component
{
    static propTypes = {
        example: PropTypes.array.isRequired,
        getExample: PropTypes.func.isRequired
    };

    componentDidMount()
    {
        this.props.getExample();
    }

    // Click Event Example
    consoleLog = (e) =>
    {
        const data = e.currentTarget.getAttribute('data-item');
        console.log(data)
    }

    render() {
        return (
            <Fragment>
                <h1>Example</h1>
                <table className="table table-striped">
                    <thead>
                        <tr>
                            <th>id</th>
                            <th>smiles</th>
                        </tr>
                    </thead>
                    <tbody>
                        {this.props.example.map(x => (
                            <tr key={x.id} data-item={x.smiles} onClick={this.consoleLog}>
                                <td>{x.id}</td>
                                <td>{x.smiles}</td>
                            </tr>
                        ))}
                    </tbody>
                </table>
            </Fragment>
        );
    };
};

const mapStateToProps = state => ({
    example: state.exampleReducer.example
});

export default connect(mapStateToProps, { getExample })(ExampleTable);