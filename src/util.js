function activate_copy_button(text_id, button_id) {
    const btn = document.getElementById(button_id);
    const txt = document.getElementById(text_id).value;
    btn.addEventListener('click', () => {
        navigator.clipboard.writeText(txt);
        btn.innerHTML = 'OK';
        setTimeout(() => (btn.innerHTML = 'Copy'), 1000);
    });
}

